#pragma once

#include "scene.h"
#include "pcg.h"
#include "transform.h"

struct Sample{
    int light_id = -1;
    PointAndNormal point;
    Sample() = default;
    Sample(int lid,PointAndNormal p):light_id(lid),point(p){};
};

struct Reservoir{
    PathVertex visible_sample;
    Sample light_sample;
    Vector3 dir_view;
    Real w_sum = Real(0);
    Real W = 0;
    int M = 0;
};


bool update(const Scene&scene,Reservoir& target,Sample y,Real w_y,pcg32_state &rng){
    target.w_sum += w_y;
    target.M++;
    if(target.w_sum > Real(0) && next_pcg32_real<Real>(rng) < (w_y / target.w_sum)){
        target.light_sample = y;
        return true;
    }
    return false;
}

bool is_hit(const Reservoir& r){
    return (r.visible_sample.material_id >= 0);
}

bool nbor_valid(const Scene&scene,const Reservoir& q,const Reservoir& n){
    if(!is_hit(q) || !is_hit(n))return false;
    Real dq = abs(xform_point(scene.camera.world_to_cam,q.visible_sample.position).z);
    Real dn = abs(xform_point(scene.camera.world_to_cam,n.visible_sample.position).z);
    // std::cout << abs(dn-dq) / dq << std::endl;
    assert(dq > 0);
    if(abs(dn-dq) / dq > 0.1){
        // std::cout << "retire pixel because of over depth" << std::endl;
        return false;
    }
    Vector3 nq = normalize(q.visible_sample.geometric_normal);
    Vector3 nn = normalize(n.visible_sample.geometric_normal);
    if(dot(nq,nn) < cos(radians(25))){
        // std::cout << "retire pixel because of over normal" << std::endl;
        return false;
    }
    // Real pdf_hat_q = 0;
    // {
    //     PointAndNormal point_on_light = q.light_sample.point;
    //     PathVertex vertex = q.visible_sample;
    //     Vector3 dir_view = q.dir_view;
    //     Vector3 dir_light = normalize(point_on_light.position - vertex.position);
    //     const Light& light = scene.lights[q.light_sample.light_id];
    //     Ray shadow_ray{vertex.position, dir_light, 
    //         get_shadow_epsilon(scene),
    //         (1 - get_shadow_epsilon(scene)) *
    //             distance(point_on_light.position, vertex.position)};
    //     Real V = occluded(scene, shadow_ray)?0.0:1.0;
    //     Real G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
    //         distance_squared(point_on_light.position, vertex.position);
    //     assert(vertex.material_id >= 0);
    //     const Material &mat = scene.materials[vertex.material_id];
    //     Spectrum f = eval(mat, dir_view, dir_light, vertex, scene.texture_pool);
    //     Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
    //     pdf_hat_q = luminance(f * L * G * V);
    // }

    // Real pdf_hat_n = 0.0;
    // {
    //     PointAndNormal point_on_light = n.light_sample.point;
    //     PathVertex vertex = q.visible_sample;
    //     Vector3 dir_view = q.dir_view;
    //     Vector3 dir_light = normalize(point_on_light.position - vertex.position);
    //     const Light& light = scene.lights[n.light_sample.light_id];
    //     Ray shadow_ray{vertex.position, dir_light, 
    //         get_shadow_epsilon(scene),
    //         (1 - get_shadow_epsilon(scene)) *
    //             distance(point_on_light.position, vertex.position)};
    //     Real V = occluded(scene, shadow_ray)?0.0:1.0;
    //     Real G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
    //         distance_squared(point_on_light.position, vertex.position);
    //     assert(vertex.material_id >= 0);
    //     const Material &mat = scene.materials[vertex.material_id];
    //     Spectrum f = eval(mat, dir_view, dir_light, vertex, scene.texture_pool);
    //     Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
    //     pdf_hat_n = luminance(f * L * G * V);
    // }
    // if(abs(pdf_hat_q-pdf_hat_n) / abs(pdf_hat_q) > 0.1)return false;
    return true;
}

Reservoir combine(const Scene& scene,const Reservoir& q,std::vector<const Reservoir*> nbors,pcg32_state &rng){
    Reservoir s = q;
    int nM = q.M;
    // bool selected = false;
    for(const auto& r:nbors){
        // calculate ph_q(ri)
        if(!nbor_valid(scene,q,*r))continue;
        PointAndNormal point_on_light = r->light_sample.point;
        PathVertex vertex = q.visible_sample;
        Vector3 dir_view = q.dir_view;
        Vector3 dir_light = normalize(point_on_light.position - vertex.position);
        assert(r->light_sample.light_id >= 0);
        const Light& light = scene.lights[r->light_sample.light_id];
        Ray shadow_ray{vertex.position, dir_light, 
            get_shadow_epsilon(scene),
            (1 - get_shadow_epsilon(scene)) *
                distance(point_on_light.position, vertex.position)};
        Real V = occluded(scene, shadow_ray)?0.0:1.0;
        Real G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
            distance_squared(point_on_light.position, vertex.position);
        assert(vertex.material_id >= 0);
        const Material &mat = scene.materials[vertex.material_id];
        Spectrum f = eval(mat, dir_view, dir_light, vertex, scene.texture_pool);
        Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
        Real pdf_hat = luminance(f * L * G);
        Real wy = pdf_hat * r->W * r->M;
        update(scene,s,r->light_sample,wy,rng);
        // selected = 
        nM += r->M;
    }
    s.M = nM;
    int Z = q.M;
    for(const auto& r:nbors){
        if(!nbor_valid(scene,q,*r))continue;
        PointAndNormal point_on_light = s.light_sample.point;
        PathVertex vertex = r->visible_sample;
        Vector3 dir_view = r->dir_view;
        Vector3 dir_light = normalize(point_on_light.position - vertex.position);
        assert(r->light_sample.light_id >= 0);
        const Light& light = scene.lights[r->light_sample.light_id];
        Ray shadow_ray{vertex.position, dir_light, 
            get_shadow_epsilon(scene),
            (1 - get_shadow_epsilon(scene)) *
                distance(point_on_light.position, vertex.position)};
        Real V = occluded(scene, shadow_ray)?0.0:1.0;
        Real G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
            distance_squared(point_on_light.position, vertex.position);
        assert(vertex.material_id >= 0);
        const Material &mat = scene.materials[vertex.material_id];
        Spectrum f = eval(mat, dir_view, dir_light, vertex, scene.texture_pool);
        Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
        Real pdf_hat = luminance(f * L * G);
        if(pdf_hat > Real(0)){
            Z += r->M;
        }
    }
    assert(Real(Z) >= 0);
    Real m = Real(1) / Real(Z);
    PointAndNormal point_on_light = q.light_sample.point;
    PathVertex vertex = q.visible_sample;
    Vector3 dir_view = q.dir_view;
    Vector3 dir_light = normalize(point_on_light.position - vertex.position);
    const Light& light = scene.lights[q.light_sample.light_id];
    Real G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
        distance_squared(point_on_light.position, vertex.position);
    Ray shadow_ray{vertex.position, dir_light, 
        get_shadow_epsilon(scene),
        (1 - get_shadow_epsilon(scene)) *
            distance(point_on_light.position, vertex.position)};
    Real V = occluded(scene, shadow_ray)?0.0:1.0;
    assert(vertex.material_id >= 0);
    const Material &mat = scene.materials[vertex.material_id];
    Spectrum f = eval(mat, dir_view, dir_light, vertex, scene.texture_pool);
    Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
    Real pdf_hat = luminance(f * L * G);
    
    s.W = Real(1) / pdf_hat * (m * s.w_sum);
    // std::cout << "combined" << std::endl;

    return s;
}

void initial_sample(const Scene& scene, Reservoir* buffer,int x,int y,pcg32_state &rng){
    int w = scene.camera.width, h = scene.camera.height;
    int buffer_idx = y * w + x;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    Vector3 dir_view = -ray.dir;
    RayDifferential ray_diff = init_ray_differential(w, h);
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);

    Reservoir rr;
    rr.dir_view = -ray.dir;
    if(vertex_.has_value()){
        PathVertex vertex = *vertex_;
        rr.visible_sample = vertex;
        if(!is_light(scene.shapes[vertex.shape_id])){
            bool is_selected = false;
            for(int i = 0;i<scene.options.initial_M || !is_selected;++i){
                // RIS
                Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
                Real light_w = next_pcg32_real<Real>(rng);
                Real shape_w = next_pcg32_real<Real>(rng);
                int light_id = sample_light(scene, light_w);
                const Light &light = scene.lights[light_id];
                PointAndNormal point_on_light =
                    sample_point_on_light(light, vertex.position, light_uv, shape_w, scene);
                
                Real l_pdf = light_pmf(scene,light_id) * pdf_point_on_light(scene.lights[light_id],point_on_light,vertex.position,scene);
                Vector3 dir_light = normalize(point_on_light.position - vertex.position);
                Ray shadow_ray{vertex.position, dir_light, 
                    get_shadow_epsilon(scene),
                    (1 - get_shadow_epsilon(scene)) *
                        distance(point_on_light.position, vertex.position)};
                Real V = occluded(scene, shadow_ray)?0.0:1.0;
                Real G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
                    distance_squared(point_on_light.position, vertex.position);
                assert(vertex.material_id >= 0);
                const Material &mat = scene.materials[vertex.material_id];
            
                Spectrum f = eval(mat, dir_view, dir_light, vertex, scene.texture_pool);
                Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
                
                Spectrum fy = f * L * G;
                Real pdf_hat = luminance(fy);
                Real w = pdf_hat / l_pdf;
                Sample ls(light_id,point_on_light);
                if(update(scene,rr,ls,w,rng)){
                    is_selected = true;
                    rr.W = Real(1) / pdf_hat * (Real(1)/Real(rr.M)) * rr.w_sum;
                }
            }
        }else{
            std::cout << "hit light" << std::endl;
            // rr.M = rr.ph = rr.w_sum = 1.0;
            // rr.f = emission(vertex,-ray.dir,scene);
            // rr.visible_sample = vertex;
            // PointAndNormal pn{vertex.position,vertex.geometric_normal};
            // rr.light_sample = {get_area_light_id(scene.shapes[vertex.shape_id]),pn};
        }
    }else{
        std::cout << "hit nothing" << std::endl;
    }
    buffer[buffer_idx] = rr;
}

Spectrum shadow_pixel(const Scene& scene, Reservoir* buffer,int x,int y,pcg32_state&rng){
    int w = scene.camera.width, h = scene.camera.height;
    int buffer_idx = y * w + x;
    Reservoir& rr = buffer[buffer_idx];
    PathVertex vertex = rr.visible_sample;
    PointAndNormal point_on_light = rr.light_sample.point;
    Vector3 dir_light = normalize(point_on_light.position - vertex.position);
    Ray shadow_ray{vertex.position, dir_light, 
        get_shadow_epsilon(scene),
        (1 - get_shadow_epsilon(scene)) *
            distance(point_on_light.position, vertex.position)};
    Real visible = 1.0;
    // debug
    if(occluded(scene,shadow_ray)){
        rr.W = Real(0);
        visible = 0.0;
    }
    return {visible,visible,visible};
}

void spatial_combine(const Scene& scene, const Reservoir* original_buffer, Reservoir* target_buffer,int x, int y,pcg32_state &rng){
    int w = scene.camera.width, h = scene.camera.height;
    int buffer_idx = y * w + x;
    // target_buffer[buffer_idx] = original_buffer[buffer_idx];
    Reservoir q = original_buffer[buffer_idx];
    Reservoir &rr = target_buffer[buffer_idx]; 
    int radius = scene.options.pixel_radius;
    int x_anchor = x - radius;
    int y_anchor = y - radius;
    int tx = 0,ty = 0;
    std::vector<const Reservoir*> nbors;
    int cnt = 0;
    while(cnt < scene.options.k){
        tx = static_cast<int>(x_anchor + next_pcg32_real<Real>(rng) * 2 * radius);
        ty = static_cast<int>(y_anchor + next_pcg32_real<Real>(rng) * 2 * radius);
        if(tx < 0 || tx >= w || ty < 0 || ty >= h)continue;
        if(!nbor_valid(scene,q,original_buffer[ty * w + tx]))continue;
        nbors.push_back(&original_buffer[ty * w + tx]);
        ++cnt;
    }
    rr = combine(scene,q,nbors,rng);
}

Spectrum shade_pixel(const Scene& scene, const Reservoir* r_buffer, int x,int y,pcg32_state &rng){
    int w = scene.camera.width, h = scene.camera.height;
    int buffer_idx = y * w + x;
    const Reservoir& rr = r_buffer[buffer_idx];
    // if(!is_hit(rr))
    PointAndNormal point_on_light = rr.light_sample.point;
    PathVertex vertex = rr.visible_sample;
    Vector3 dir_light = normalize(point_on_light.position - vertex.position);
    Ray shadow_ray{vertex.position, dir_light, 
        get_shadow_epsilon(scene),
        (1 - get_shadow_epsilon(scene)) *
            distance(point_on_light.position, vertex.position)};
    Real V = occluded(scene, shadow_ray)?0.0:1.0;
    Real G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
        distance_squared(point_on_light.position, vertex.position);
    assert(vertex.material_id >= 0);
    const Material &mat = scene.materials[vertex.material_id];

    Spectrum f = eval(mat, rr.dir_view, dir_light, vertex, scene.texture_pool);
    assert(rr.light_sample.light_id >= 0);
    Spectrum L = emission(scene.lights[rr.light_sample.light_id], -dir_light, Real(0), point_on_light, scene);
    
    Spectrum fy = f * L * G * V;
    return fy * rr.W;
}