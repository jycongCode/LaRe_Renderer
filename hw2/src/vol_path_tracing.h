#pragma once

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = {Real(0),Real(0)};
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    if(vertex_){
        PathVertex vertex = *vertex_;
        Real t_hit = distance(vertex.position,ray.org);
        int exterior_medium_id = get_exterior_medium_id(scene.shapes[vertex.shape_id]);
        Medium medium = scene.media[exterior_medium_id];
        Spectrum sigma_a = get_sigma_a(medium,vertex.position);
        Spectrum transmittance = exp(-sigma_a * t_hit);
        Spectrum Le(Real(0),Real(0),Real(0));
        if(is_light(scene.shapes[vertex.shape_id])){
            Le = emission(vertex, -ray.dir, scene);
        }
        return transmittance * Le;
    }
    return make_zero_spectrum();
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = {Real(0),Real(0)};
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    bool is_hit = vertex_.has_value();
    PathVertex vertex;
    Real t_hit = ray.tfar;
    Medium medium;
    if(is_hit){
        vertex = *vertex_;
        t_hit = distance(vertex.position,ray.org);
        int exterior_medium_id = get_exterior_medium_id(scene.shapes[vertex.shape_id]);
        medium = scene.media[exterior_medium_id];
    }else{
        int exterior_medium_id = scene.camera.medium_id;
        medium = scene.media[exterior_medium_id];
    }
    Spectrum sigma_a = get_sigma_a(medium,vertex.position);
    Spectrum sigma_s = get_sigma_s(medium,vertex.position);
    Spectrum sigma_t = sigma_a + sigma_s;
    Real u = next_pcg32_real<Real>(rng);
    Real t = -log(1-u)/sigma_t.x;
    if (!is_hit || t < t_hit){
        Spectrum trans_pdf = exp(-sigma_t * t)*sigma_t;
        Spectrum transmittance = exp(-sigma_t * t);
        Vector3 p = ray.org + t * ray.dir;
        PhaseFunction phasefunc = get_phase_function(medium);
        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);
        int light_id = sample_light(scene, light_w);
        const Light &light = scene.lights[light_id];
        PointAndNormal point_on_light =
            sample_point_on_light(light, p, light_uv, shape_w, scene);
        Vector3 dir_light = normalize(point_on_light.position-p);
        Spectrum rho = eval(phasefunc,-ray.dir,dir_light);
        Spectrum light_emission = emission(light, point_on_light.normal, Real(0), point_on_light, scene);
        Ray shadow_ray{p, dir_light, 
                            0,
                            distance(p,point_on_light.position)};
        Real G = 1.0;
        if(occluded(scene,shadow_ray)){
            G = 0.0;
        }
        Spectrum L_s1_estimate = rho * light_emission * exp(-sigma_t * distance(p,point_on_light.position)) \
                        * dot(-dir_light,point_on_light.normal) / distance_squared(p,point_on_light.position) * G;
        Real L_s1_pdf = light_pmf(scene,light_id) * pdf_point_on_light(light,point_on_light,p,scene);
        return (transmittance / trans_pdf) * sigma_s * (L_s1_estimate / L_s1_pdf); 
    }else{
        Spectrum trans_pdf = exp(-sigma_t * t_hit);
        Spectrum transmittance = exp(-sigma_t * t_hit);
        Spectrum Le(Real(0),Real(0),Real(0));
        if(is_light(scene.shapes[vertex.shape_id])){
            Le = emission(vertex,-ray.dir,scene);
        }
        return (transmittance / trans_pdf) * Le;
    }
    return make_zero_spectrum();
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
int update_medium(PathVertex vertex,Ray ray,int& medium){
    if(vertex.exterior_medium_id != vertex.interior_medium_id)
        medium = dot(ray.dir,vertex.geometric_normal)>Real(0) ? vertex.exterior_medium_id : vertex.interior_medium_id;
    return medium;
}

Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = {Real(0),Real(0)};

    int current_medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput(1.0,1.0,1.0);
    Spectrum radiance(0.0,0.0,0.0);
    int bounces = 0;
    int max_depth = scene.options.max_depth;
    while(true){
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene,ray,ray_diff);
        PathVertex vertex;
        if(vertex_.has_value()){
            vertex = *vertex_;
        }
        Spectrum transmittance(1.0,1.0,1.0);
        Spectrum trans_pdf(1.0,1.0,1.0);

        if(current_medium_id!=-1){
            Vector3 p = ray.org;
            if(vertex_.has_value()){
                p = vertex.position;
            }
            Medium current_medium = scene.media[current_medium_id];
            Spectrum sigma_t = get_sigma_a(current_medium,p) + get_sigma_s(current_medium,p);
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1-u)/sigma_t.x;
            transmittance = exp(-sigma_t * t);
            trans_pdf = sigma_t * exp(-sigma_t * t);
            scatter = true;
            Real t_hit = distance(vertex.position,ray.org);
            if(vertex_.has_value()){
                if(t > t_hit){
                    scatter = false;
                    transmittance = exp(-sigma_t * t_hit);
                    trans_pdf = exp(-sigma_t * t_hit);
                }
            }
            if(vertex_.has_value()&&vertex.material_id == -1&&t>t_hit){
                ray.org += t_hit * ray.dir;
            }else{
                ray.org += t * ray.dir;
            }            
        }
        current_path_throughput *= (transmittance/trans_pdf);

        if(!scatter){
            Spectrum Le(0.0,0.0,0.0);
            if(vertex_.has_value() && is_light(scene.shapes[(*vertex_).shape_id])){
                Le = emission((*vertex_),-ray.dir,scene);
            }
            radiance += current_path_throughput * Le;
        }

        if(bounces == max_depth-1 && max_depth != -1){
            break;
        }

        if(!scatter && vertex_.has_value()){
            if(vertex.material_id == -1){
                update_medium(vertex,ray,current_medium_id);
                ++bounces;
                continue;
            }
        }

        if(scatter){
            PhaseFunction phasefunc = get_phase_function(scene.media[current_medium_id]);
            Vector2 rnd_uv(next_pcg32_real<Real>(rng),next_pcg32_real<Real>(rng));
            std::optional<Vector3> new_dir_ = sample_phase_function(phasefunc,-ray.dir,rnd_uv);
            Vector3 new_dir = ray.dir;
            if(new_dir_.has_value()){
                new_dir = *new_dir_;
            }
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id],ray.org);
            current_path_throughput *= (eval(phasefunc,-ray.dir,new_dir)/pdf_sample_phase(phasefunc,-ray.dir,new_dir))*sigma_s;
            ray.dir = new_dir;
        }else{
            break;
        }

        Real rr_prob = 1.0;
        if(bounces >= scene.options.rr_depth){
            rr_prob = min(current_path_throughput.x,0.95);
            if(next_pcg32_real<Real>(rng) > rr_prob){
                break;
            }else{
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
    }
    return radiance;

}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase get_shadow_epsilon(scene)function sampling
// still no surface lighting

Spectrum next_event_estimation(bool is_volume,PathVertex vertex, Ray ray, int medium,int bounces, pcg32_state &rng, const Scene&scene){
    Vector3 p = vertex.position;
    Vector3 start_p = vertex.position;
    // sample point on light
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];
    PointAndNormal light_isect =
        sample_point_on_light(light, p, light_uv, shape_w, scene);
    Vector3 p_prime = light_isect.position;
    
    // params
    Spectrum T_light = make_const_spectrum(1);
    Spectrum p_trans_dir = make_const_spectrum(1);
    Spectrum p_trans_nee = make_const_spectrum(1);
    int nee_bounces = 0;
    int current_medium_id = medium;
    int max_depth = scene.options.max_depth;
    Vector3 dir_light = normalize(p_prime-p);
    Vector3 dir_view = -ray.dir;
    RayDifferential ray_diff{0.0,0.0};
    while(true){
        Ray nee_ray{p, dir_light, 
                    get_shadow_epsilon(scene),
                    (1 - get_shadow_epsilon(scene)) *
                        distance(p_prime, p)};

        std::optional<PathVertex> shadow_vertex_ = intersect(scene,nee_ray,ray_diff);
        Real next_t = distance(p_prime,p);
        PathVertex shadow_vertex;
        bool is_hit = false;
        if(shadow_vertex_.has_value()){
            is_hit = true;
            shadow_vertex = *shadow_vertex_;
            next_t = distance(shadow_vertex.position,p);
        }

        if(current_medium_id != -1){
            Real u = next_pcg32_real<Real>(rng);
            int channel = std::clamp(static_cast<int>(u*3),0,2);
            int iteration = 0;
            Real accum_t = 0.0;
            Spectrum majorant = get_majorant(scene.media[current_medium_id],nee_ray);
            while(true){
                if(majorant[channel] <= 0)break;
                if(iteration >= scene.options.max_null_collisions)break;
                Real t = -log(1.0-next_pcg32_real<Real>(rng)) / majorant[channel];
                Real dt = next_t - accum_t;
                accum_t = min(accum_t + t,next_t);
                Spectrum sigma_a = get_sigma_a(scene.media[current_medium_id],nee_ray.org+accum_t * nee_ray.dir);
                Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id],nee_ray.org+accum_t * nee_ray.dir);
                Spectrum sigma_t = sigma_a + sigma_s;
                Spectrum sigma_n = majorant - sigma_t;
                if(t < dt){
                    T_light *= exp(-majorant * t) * sigma_n / max(majorant);
                    p_trans_nee *= exp(-majorant * t) * majorant / max(majorant);
                    Spectrum real_prob = sigma_t / majorant;
                    p_trans_dir *= exp(-majorant * t) * majorant * (1 - real_prob[channel]) / max(majorant);
                    if(max(T_light) <= 0)break;
                }else{
                    T_light *= exp(-majorant * dt);
                    p_trans_nee *= exp(-majorant * dt);
                    p_trans_dir *= exp(-majorant * dt);
                    break;
                }
                iteration += 1;
            }
        }

        if(!is_hit)break;
        else{
            if(shadow_vertex.material_id >= 0){
                return Spectrum(0.0,0.0,0.0);
            }
            nee_bounces += 1;
            if(max_depth != -1 && bounces + nee_bounces + 1 >= max_depth){
                return Spectrum(0.0,0.0,0.0);
                // break;
            }
            current_medium_id = update_medium(shadow_vertex,nee_ray,current_medium_id);
            p += next_t * nee_ray.dir;
        }
    }

    if(max(T_light) > Real(0)){
        if(is_volume){
            PhaseFunction phasefunc = get_phase_function(scene.media[medium]);
            Spectrum sigma_s = get_sigma_s(scene.media[medium],start_p);
            Spectrum rho = eval(phasefunc,dir_view,dir_light) * sigma_s;
            Spectrum Le = emission(light,-dir_light,0.0,light_isect,scene);
            Real G = fabs(dot(-dir_light,light_isect.normal)) / distance_squared(start_p,p_prime);
            
            Real pdf_nee = average(p_trans_nee) * light_pmf(scene,light_id) * pdf_point_on_light(light,light_isect,start_p,scene);
            Spectrum contrib = T_light * G * rho * Le / pdf_nee;
            Real pdf_phase = pdf_sample_phase(phasefunc,dir_view,dir_light) * G * average(p_trans_dir);
            Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
            return w * contrib;
        }else{
            const Material& mat = scene.materials[vertex.material_id];
            Spectrum rho = eval(mat, dir_view, dir_light, vertex, scene.texture_pool);
            Spectrum Le = emission(light,-dir_light,Real(0),light_isect,scene);
            Real G = fabs(dot(-dir_light,light_isect.normal)) / distance_squared(start_p,p_prime);
            
            Real pdf_nee = average(p_trans_nee) * light_pmf(scene,light_id) * pdf_point_on_light(light,light_isect,start_p,scene);
            Spectrum contrib = T_light * G * rho * Le / pdf_nee;
            Real pdf_bsdf = pdf_sample_bsdf(mat,dir_view,dir_light,vertex,scene.texture_pool) * G * average(p_trans_dir);
            Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_bsdf * pdf_bsdf);
            return w * contrib;
        }
    }
    return {0.0,0.0,0.0};
}

Spectrum volume_nee(const Scene &scene,Vector3 p, int current_medium_id,Vector3 ray_dir,int bounces,pcg32_state &rng){
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];
    PointAndNormal light_isect =
        sample_point_on_light(light, p, light_uv, shape_w, scene);
    Vector3 p_prime = light_isect.position;
    Spectrum T_light(1.0,1.0,1.0);
    int shadow_medium_id = current_medium_id;
    int shadow_bounces = 0;
    Vector3 p_trans_dir{1.0,1.0,1.0};
    int max_depth = scene.options.max_depth;
    Vector3 dir_light = normalize(p_prime - p);
    Vector3 start_p = p;
    while(true){
        Ray shadow_ray{p, normalize(p_prime-p), 
                               get_shadow_epsilon(scene),
                               (1 - get_shadow_epsilon(scene)) *
                                   distance(p_prime, p)};
        RayDifferential shadow_ray_diff{0.0,0.0};
        std::optional<PathVertex> vertex_ = intersect(scene,shadow_ray,shadow_ray_diff);
        PathVertex vertex;
        bool is_hit = false;
        Real next_t = distance(p_prime,p);
        if(vertex_.has_value()){
            is_hit = true;
            vertex = *vertex_;
            next_t = distance(vertex.position,shadow_ray.org);
        }
        
        if(shadow_medium_id != -1){
            Spectrum sigma_a = get_sigma_a(scene.media[shadow_medium_id],p);
            Spectrum sigma_s = get_sigma_s(scene.media[shadow_medium_id],p);
            Spectrum sigma_t = sigma_a + sigma_s;
            T_light *= exp(-sigma_t * next_t);
            p_trans_dir *= exp(-sigma_t.x * next_t);
        }

        if(!is_hit){
            break;
        }else{
            if(vertex.material_id >= 0){
                return Spectrum(0.0,0.0,0.0);
            }
            shadow_bounces += 1;
            if(max_depth != -1 && bounces + shadow_bounces + 1 >= max_depth){
                // return Spectrum(0.0,0.0,0.0);
                break;
            }
            shadow_medium_id = update_medium(vertex,shadow_ray,shadow_medium_id);
            p += next_t * shadow_ray.dir;
        }
    }
    if(T_light.x > 0.0 && T_light.y > 0.0 && T_light.z > 0.0){
        PhaseFunction phasefunc = get_phase_function(scene.media[current_medium_id]);
        Spectrum rho = eval(phasefunc,-ray_dir,dir_light);
        Spectrum Le = emission(scene.lights[light_id],-dir_light,1.0,light_isect,scene);
        Real pdf_nee = light_pmf(scene,light_id) * pdf_point_on_light(scene.lights[light_id],light_isect,start_p,scene);
        Real G = fabs(dot(-dir_light,light_isect.normal)) / distance_squared(start_p,p_prime);
        Spectrum contrib = T_light * G * rho * Le / pdf_nee;
        Vector3 pdf_phase = pdf_sample_phase(phasefunc,-ray_dir,dir_light) * G * p_trans_dir;
        Vector3 w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
        return w * contrib;
    }
    return {0.0,0.0,0.0};
}

Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = {Real(0),Real(0)};

    int current_medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput(1.0,1.0,1.0);
    Spectrum radiance(0.0,0.0,0.0);
    int bounces = 0;
    Real dir_pdf = 0.0;
    bool volume_never_scatter = true;
    Spectrum multi_trans_pdf(1.0,1.0,1.0);
    Vector3 nee_p_cache;
   
    int max_depth = scene.options.max_depth;
    while(true){
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene,ray,ray_diff);
        PathVertex vertex;
        bool is_hit = false;
        if(vertex_.has_value()){
            vertex = *vertex_;
            is_hit = true;
        }
        Spectrum transmittance(1.0,1.0,1.0);
        Spectrum trans_pdf(1.0,1.0,1.0);
        Vector3 p = ray.org;
        // not travelling in vaccum
        if(current_medium_id!=-1){
            Medium current_medium = scene.media[current_medium_id];
            
            Spectrum sigma_t = get_sigma_a(current_medium,p) + get_sigma_s(current_medium,p);
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1-u)/sigma_t.x;
            transmittance = exp(-sigma_t * t);
            trans_pdf = sigma_t * exp(-sigma_t * t);
            scatter = true;

            Real t_hit = 0.0;
            if(is_hit){
                t_hit = distance(vertex.position,ray.org);
                if(t > t_hit){
                    scatter = false;
                    transmittance = exp(-sigma_t * t_hit);
                    trans_pdf = exp(-sigma_t * t_hit);
                    if(vertex.material_id == -1){
                        ray.org += t_hit * ray.dir;
                    }
                }else{
                    scatter = true;
                    ray.org += t * ray.dir;
                }
            }else{
                scatter = true;
                transmittance = exp(-sigma_t * t);
                trans_pdf = sigma_t * exp(-sigma_t * t);
                ray.org += t * ray.dir;
            }           
        }else{
            if(is_hit){
                ray.org += distance(ray.org,vertex.position) * ray.dir;
            }
        }
        current_path_throughput *= (transmittance/trans_pdf);
        if(scatter){
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id],ray.org);
            radiance += current_path_throughput * volume_nee(scene,ray.org,current_medium_id,ray.dir,bounces,rng) * sigma_s;
            // update latest nee position
            nee_p_cache = p;
            multi_trans_pdf = {1.0,1.0,1.0};
        }

        if(!scatter && is_hit){
            multi_trans_pdf *= trans_pdf;
            // hit light or hit index-match
            if(volume_never_scatter){
                Spectrum Le(0.0,0.0,0.0);
                if(is_light(scene.shapes[vertex.shape_id])){
                    Le = emission(vertex,-ray.dir,scene);
                }
                radiance += current_path_throughput * Le;
            }else{
                // hit light
                if(is_light(scene.shapes[vertex.shape_id])){
                    Vector3 light_point = vertex.position;
                    int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    PointAndNormal light_isect{vertex.position,vertex.geometric_normal};
                    Real pdf_nee = light_pmf(scene,light_id) * pdf_point_on_light(scene.lights[light_id],light_isect,nee_p_cache,scene);
                    Real G = fabs(dot(-ray.dir,light_isect.normal)) / distance_squared(nee_p_cache,vertex.position);
                    Vector3 dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                    Vector3 w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    radiance += current_path_throughput * emission(vertex,-ray.dir,scene) * w;
                }
            }
        }

        if(bounces == max_depth-1 && max_depth != -1){
            break;
        }

        if(!scatter && is_hit){
            if(vertex.material_id == -1){
                update_medium(vertex,ray,current_medium_id);
                bounces += 1;
                continue;
            }
        }

        if(scatter){
            volume_never_scatter = false;
            PhaseFunction phasefunc = get_phase_function(scene.media[current_medium_id]);
            Vector2 rnd_uv(next_pcg32_real<Real>(rng),next_pcg32_real<Real>(rng));
            std::optional<Vector3> new_dir_ = sample_phase_function(phasefunc,-ray.dir,rnd_uv);
            if(!new_dir_)break;
            Vector3 new_dir = *new_dir_;
            dir_pdf = pdf_sample_phase(phasefunc,-ray.dir,new_dir);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id],ray.org);
            current_path_throughput *= (eval(phasefunc,-ray.dir,new_dir)/dir_pdf)*sigma_s;
            ray.dir = new_dir;
        }else{
            break;
        }

        Real rr_prob = 1.0;
        if(bounces >= scene.options.rr_depth){
            rr_prob = min(current_path_throughput.x,0.95);
            if(next_pcg32_real<Real>(rng) > rr_prob){
                break;
            }else{
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
    }
    return radiance;
}


Spectrum surface_nee(const Scene &scene,PathVertex vertex, int current_medium_id,Vector3 ray_dir,int bounces,pcg32_state &rng){
    Vector3 p = vertex.position;
    const Material& mat = scene.materials[vertex.material_id];
    
    // sample point on light
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];
    PointAndNormal light_isect =
        sample_point_on_light(light, p, light_uv, shape_w, scene);
    Vector3 p_prime = light_isect.position;
    
    // params
    Spectrum T_light(1.0,1.0,1.0);
    int shadow_medium_id = current_medium_id;
    int shadow_bounces = 0;
    Vector3 p_trans_dir{1.0,1.0,1.0};
    int max_depth = scene.options.max_depth;
    // Ray shadow_ray{p,normalize(p_prime-p),0,distance(p,p_prime)-0.01};
    Vector3 dir_light = normalize(p_prime-p);
    Vector3 start_p = vertex.position;
    while(true){
        Ray shadow_ray{p, normalize(p_prime-p), 
                               get_shadow_epsilon(scene),
                               (1 - get_shadow_epsilon(scene)) *
                                   distance(p_prime, p)};
        RayDifferential shadow_ray_diff{0.0,0.0};
        std::optional<PathVertex> shadow_vertex_ = intersect(scene,shadow_ray,shadow_ray_diff);
        Real next_t = distance(p_prime,p);
        PathVertex shadow_vertex;
        bool is_hit = false;
        if(shadow_vertex_.has_value()){
            is_hit = true;
            shadow_vertex = *shadow_vertex_;
            next_t = distance(vertex.position,p);
        }

        if(shadow_medium_id != -1){
            Spectrum sigma_a = get_sigma_a(scene.media[shadow_medium_id],p);
            Spectrum sigma_s = get_sigma_s(scene.media[shadow_medium_id],p);
            Spectrum sigma_t = sigma_a + sigma_s;
            T_light *= exp(-sigma_t * next_t);
            p_trans_dir *= exp(-sigma_t.x * next_t);
        }

        if(!is_hit){
            break;
        }else{
            if(shadow_vertex.material_id >= 0){
                return Spectrum(0.0,0.0,0.0);
            }
            shadow_bounces += 1;
            if(max_depth != -1 && bounces + shadow_bounces + 1 >= max_depth){
                return Spectrum(0.0,0.0,0.0);
            }

            shadow_medium_id = update_medium(shadow_vertex,shadow_ray,shadow_medium_id);
            p = p + next_t * shadow_ray.dir;
        }
    }
    if(T_light.x > 0.0 && T_light.y > 0.0 && T_light.z > 0.0){
        // PhaseFunction phasefunc = get_phase_function(scene.media[current_medium_id]);
        Spectrum rho = eval(mat, -ray_dir, dir_light, vertex, scene.texture_pool);
        // Spectrum rho = eval(phasefunc,-ray_dir,shadow_ray.dir);
        Spectrum Le = emission(scene.lights[light_id],-dir_light,1.0,light_isect,scene);
        Real pdf_nee = light_pmf(scene,light_id) * pdf_point_on_light(scene.lights[light_id],light_isect,start_p,scene);
        Real G = fabs(dot(-dir_light,light_isect.normal)) / distance_squared(start_p,p_prime);
        Spectrum contrib = T_light * G * rho * Le / pdf_nee;
    
        Vector3 pdf_phase = pdf_sample_bsdf(mat,-ray_dir,dir_light,vertex,scene.texture_pool) * G * p_trans_dir;
        Vector3 w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
        return w * contrib;
    }
    return {0.0,0.0,0.0};
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff{Real(0),Real(0)};

    int current_medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput(1.0,1.0,1.0);
    Spectrum radiance(0.0,0.0,0.0);
    int bounces = 0;
    Real dir_pdf = 0.0;
    Spectrum multi_trans_pdf(1.0,1.0,1.0);
    Vector3 nee_p_cache;
    int max_depth = scene.options.max_depth;
    bool never_scatter = true;
    
    while(true){
        bool scatter = false;
        bool is_surface = false;
        std::optional<PathVertex> vertex_ = intersect(scene,ray,ray_diff);
        PathVertex vertex;
        bool is_hit = false;
        if(vertex_.has_value()){
            vertex = *vertex_;
            is_hit = true;
        }
        Spectrum transmittance(1.0,1.0,1.0);
        Spectrum trans_pdf(1.0,1.0,1.0);
        Vector3 p = ray.org;
        // not travelling in vaccum
        if(current_medium_id!=-1){
            const Medium& current_medium = scene.media[current_medium_id];
            Spectrum sigma_t = get_sigma_a(current_medium,p) + get_sigma_s(current_medium,p);
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1-u)/sigma_t.x;
            transmittance = exp(-sigma_t * t);
            trans_pdf = sigma_t * exp(-sigma_t * t);
            Real t_hit = infinity<Real>();
            if(is_hit){
                t_hit = distance(vertex.position,p);
            }
            if(t > t_hit){
                if(vertex.material_id >= 0 && !is_light(scene.shapes[vertex.shape_id])){
                    scatter = true;
                    is_surface = true;
                }
                transmittance = exp(-sigma_t * t_hit);
                trans_pdf = exp(-sigma_t * t_hit);
                ray.org += t_hit * ray.dir;
            }else{
                scatter = true;
                is_surface = false;
                transmittance = exp(-sigma_t * t);
                trans_pdf = sigma_t * exp(-sigma_t * t);
                ray.org += t * ray.dir;
            }           
        }else{
            if(is_hit){
                if(vertex.material_id >= 0){
                    if(!is_light(scene.shapes[vertex.shape_id])){    
                        scatter = true;
                        is_surface = true;
                    }
                }
                ray.org += distance(ray.org,vertex.position) * ray.dir;
            }
        }
        ray = Ray{ray.org,ray.dir,get_intersection_epsilon(scene),infinity<Real>()};
        current_path_throughput *= (transmittance/trans_pdf);
        multi_trans_pdf *= trans_pdf;
        if(!scatter && is_hit){
            if(is_light(scene.shapes[vertex.shape_id])){
                if(never_scatter){
                    Spectrum Le = emission(vertex,-ray.dir,scene);
                    radiance += current_path_throughput * Le;
                }else{
                    Vector3 light_point = vertex.position;
                    int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    PointAndNormal light_isect{vertex.position,vertex.geometric_normal};
                    Real pdf_nee = light_pmf(scene,light_id) * pdf_point_on_light(scene.lights[light_id],light_isect,nee_p_cache,scene);
                    Real G = fmax(dot(-ray.dir,light_isect.normal),Real(0)) / distance_squared(nee_p_cache,vertex.position);
                    Vector3 dir_pdf_ = dir_pdf * multi_trans_pdf * G; 
                    Vector3 w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    radiance += current_path_throughput * emission(vertex,-ray.dir,scene) * w;
                }
            }
        }

        if(bounces == max_depth-1 && max_depth != -1){
            break;
        }

        if(!scatter && is_hit){
            if(vertex.material_id == -1){
                current_medium_id = update_medium(vertex,ray,current_medium_id);
                bounces += 1;
                continue;
            }
        }

        if(scatter){
            never_scatter = false;
            multi_trans_pdf = {1.0,1.0,1.0}; 
            if(is_surface){
                Vector3 p = vertex.position;
                nee_p_cache = ray.org;
                // sample bsdf
                const Material& mat = scene.materials[vertex.material_id];
                Vector2 rnd_uv(next_pcg32_real<Real>(rng),next_pcg32_real<Real>(rng));
                Real rnd_w = next_pcg32_real<Real>(rng);
                std::optional<BSDFSampleRecord> bsdf_sample_ =
                sample_bsdf(mat,
                            -ray.dir,
                            vertex,
                            scene.texture_pool,
                            rnd_uv,
                            rnd_w);
                if(!bsdf_sample_){
                    // std::cout << "here" << std::endl;
                    break;
                }
                const BSDFSampleRecord& bsdf_sample = *bsdf_sample_;
                // radiance += current_path_throughput * next_event_estimation(false,vertex,ray,current_medium_id,bounces,rng,scene);
                // reflect or refract
                if (bsdf_sample.eta == 0) {
                    ray_diff.spread = reflect(ray_diff, vertex.mean_curvature, bsdf_sample.roughness);
                } else {
                    ray_diff.spread = refract(ray_diff, vertex.mean_curvature, bsdf_sample.eta, bsdf_sample.roughness);
                    current_medium_id = update_medium(vertex,ray,current_medium_id);
                }
                // scatter
                Vector3 new_dir = bsdf_sample.dir_out;
                dir_pdf = pdf_sample_bsdf(mat,-ray.dir,new_dir,vertex,scene.texture_pool);
                current_path_throughput *= (eval(mat,-ray.dir,new_dir,vertex,scene.texture_pool)/dir_pdf);
                ray.dir = new_dir;
                ray = Ray{ray.org,new_dir,get_intersection_epsilon(scene),infinity<Real>()};
            }else{
                nee_p_cache = ray.org;
                Vector3 p = ray.org;
                PathVertex scatter_vertex;
                scatter_vertex.position = p;
                PhaseFunction phasefunc = get_phase_function(scene.media[current_medium_id]);
                Vector2 rnd_uv(next_pcg32_real<Real>(rng),next_pcg32_real<Real>(rng));
                std::optional<Vector3> new_dir_ = sample_phase_function(phasefunc,-ray.dir,rnd_uv);
                if(!new_dir_){
                    break;
                }
                Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id],p);
                // radiance += current_path_throughput * next_event_estimation(true,scatter_vertex,ray,current_medium_id,bounces,rng,scene) * sigma_s;
                // radiance += current_path_throughput * volume_nee(scene,p,current_medium_id,ray.dir,bounces,rng);
                Vector3 new_dir = *new_dir_;
                dir_pdf = pdf_sample_phase(phasefunc,-ray.dir,new_dir);
                current_path_throughput *= (eval(phasefunc,-ray.dir,new_dir)/dir_pdf)*sigma_s;
                ray.dir = new_dir;
                ray = Ray{ray.org,new_dir,get_intersection_epsilon(scene),infinity<Real>()};
            }
        }else{
            break;
        }

        Real rr_prob = 1.0;
        if(bounces >= scene.options.rr_depth){
            rr_prob = min(current_path_throughput.x,0.95);
            if(next_pcg32_real<Real>(rng) > rr_prob){
                break;
            }else{
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
    }
    return radiance;
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    // std::cout << "here" << std::endl;
    
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff{Real(0),Real(0)};

    int current_medium_id = scene.camera.medium_id;

    Spectrum current_path_throughput(1.0,1.0,1.0);
    Spectrum radiance(0.0,0.0,0.0);
    int bounces = 0;
    Real dir_pdf = 0.0;
    Spectrum multi_trans_pdf = make_const_spectrum(1);
    Spectrum multi_trans_nee = make_const_spectrum(1);
    Vector3 nee_p_cache;
    int max_depth = scene.options.max_depth;
    bool never_scatter = true;
    while(true){
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene,ray,ray_diff);
        PathVertex vertex;
        bool is_hit = false;
        Real t_hit = infinity<Real>();
        if(vertex_.has_value()){
            vertex = *vertex_;
            is_hit = true;
            t_hit = distance(vertex.position,ray.org);
        }
        Spectrum transmittance = make_const_spectrum(1);
        Spectrum trans_dir_pdf = make_const_spectrum(1);
        Spectrum trans_nee_pdf = make_const_spectrum(1);
        // not travelling in vaccum
        if(current_medium_id >= 0){
            // if(current_medium_id < scene.media.size())
            Spectrum majorant = get_majorant(scene.media[current_medium_id],ray);

            // sample a channel for sampling
            Real u = next_pcg32_real<Real>(rng);
            int channel = std::clamp(int(u*3),0,2);
            Real accum_t = 0.0;
            int iteration = 0;
            while(true){ 
                if(majorant[channel] <= Real(0))break;
                if(iteration >= scene.options.max_null_collisions)break;
                Real t = -log(1.0-next_pcg32_real<Real>(rng)) / majorant[channel];
                Real dt = t_hit - accum_t;
                accum_t = min(accum_t + t,t_hit);
                if(t < dt){
                    Vector3 p = ray.org + accum_t * ray.dir;
                    Spectrum sigma_a = get_sigma_a(scene.media[current_medium_id],p);
                    Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id],p);
                    Spectrum sigma_t = sigma_a + sigma_s;
                    Spectrum sigma_n = majorant - sigma_t;
                    // sample from real/fake particle events
                    Spectrum real_prob = sigma_t / majorant;
                    if(next_pcg32_real<Real>(rng) < real_prob[channel]){
                        // hit "real" particle
                        scatter = true;
                        never_scatter = false;
                        transmittance *= exp(-majorant * t) / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * real_prob / max(majorant);
                        vertex.position = ray.org + accum_t * ray.dir;
                        vertex.exterior_medium_id = current_medium_id;
                        vertex.interior_medium_id = current_medium_id;
                        break;
                    }else{
                        // hit a "fake" particle;
                        transmittance *= exp(-majorant * t) * sigma_n / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * (1-real_prob) / max(majorant);
                        trans_nee_pdf *= exp(-majorant * t) * majorant / max(majorant);
                    }
                }else{
                    transmittance *= exp(-majorant * dt);
                    trans_dir_pdf *= exp(-majorant * dt);
                    trans_nee_pdf *= exp(-majorant * dt);
                    break;
                }
                iteration += 1;
            }         
        }
        current_path_throughput *= (transmittance/average(trans_dir_pdf));
        multi_trans_pdf *= trans_dir_pdf;
        multi_trans_nee *= trans_nee_pdf;
        if(!scatter && !is_hit)break;
        if(!scatter && is_hit){
            if(is_light(scene.shapes[vertex_->shape_id])){
                if(never_scatter){
                    Spectrum Le = emission(*vertex_,-ray.dir,scene);
                    radiance += current_path_throughput * Le;
                    break;
                }else{
                    Vector3 light_point = vertex.position;
                    int light_id = get_area_light_id(scene.shapes[vertex_->shape_id]);
                    PointAndNormal light_isect{vertex_->position,vertex_->geometric_normal};
                    Vector3 dir_light = normalize(vertex_->position-nee_p_cache);
                    Real pdf_nee = average(multi_trans_nee) * light_pmf(scene,light_id) * pdf_point_on_light(scene.lights[light_id],light_isect,nee_p_cache,scene);
                    Real G = fabs(dot(dir_light,light_isect.normal)) / distance_squared(nee_p_cache,light_isect.position);
                    Real dir_pdf_ = average(multi_trans_pdf) * dir_pdf * G; 
                    Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    radiance += current_path_throughput * emission(vertex,-ray.dir,scene) * w;
                    break;
                }
            }
        }
        
        if(bounces >= max_depth-1 && max_depth != -1){
            break;
        }

        if(!scatter && is_hit){
            if(vertex.material_id == -1){
                ray = Ray{vertex.position,ray.dir,get_intersection_epsilon(scene),infinity<Real>()};
                current_medium_id = update_medium(*vertex_,ray,current_medium_id);
                bounces += 1;
                continue;
            }
        }

        never_scatter = false;
        nee_p_cache = vertex.position;
        multi_trans_pdf = make_const_spectrum(1);
        multi_trans_nee = make_const_spectrum(1);
        Vector3 dir_in = -ray.dir;
        radiance += current_path_throughput * next_event_estimation(scatter,vertex,ray,current_medium_id,bounces,rng,scene);

        if(scatter){
            // nee
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id],vertex.position);
            // sample new direction
            PhaseFunction phasefunc = get_phase_function(scene.media[current_medium_id]);
            Vector2 rnd_uv(next_pcg32_real<Real>(rng),next_pcg32_real<Real>(rng));
            std::optional<Vector3> new_dir_ = sample_phase_function(phasefunc,-ray.dir,rnd_uv);
            if(!new_dir_){
                break;
            }
            Vector3 new_dir = *new_dir_;
            // update current_path_throughput
            dir_pdf = pdf_sample_phase(phasefunc,-ray.dir,new_dir);
            current_path_throughput *= sigma_s * (eval(phasefunc,-ray.dir,new_dir)/dir_pdf);
            ray = Ray{vertex.position,new_dir,get_intersection_epsilon(scene),infinity<Real>()};
        }else{
            // sample bsdf
            const Material& mat = scene.materials[vertex.material_id];
            Vector2 rnd_uv(next_pcg32_real<Real>(rng),next_pcg32_real<Real>(rng));
            Real rnd_w = next_pcg32_real<Real>(rng);
            std::optional<BSDFSampleRecord> bsdf_sample_ =
            sample_bsdf(mat,
                        -ray.dir,
                        vertex,
                        scene.texture_pool,
                        rnd_uv,
                        rnd_w);
            if(!bsdf_sample_){
                break;
            }
            const BSDFSampleRecord& bsdf_sample = *bsdf_sample_;
            Vector3 new_dir = bsdf_sample.dir_out;
            dir_pdf = pdf_sample_bsdf(mat,-ray.dir,new_dir,vertex,scene.texture_pool);
            current_path_throughput *= (eval(mat,-ray.dir,new_dir,vertex,scene.texture_pool)/dir_pdf);
            ray = Ray{vertex.position,new_dir,get_intersection_epsilon(scene),infinity<Real>()};
        }
        current_medium_id = update_medium(vertex,ray,current_medium_id);
        
        Real rr_prob = 1.0;
        if(bounces >= scene.options.rr_depth){
            rr_prob = min(max(current_path_throughput),0.95);
            if(next_pcg32_real<Real>(rng) > rr_prob){
                break;
            }else{
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
    }
    return radiance;
}
