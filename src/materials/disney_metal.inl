#include "../microfacet.h"


Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real anistropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum baseColor = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3 half = normalize(dir_in + dir_out);
    Real h_wout = fmax(dot(half, dir_out), Real(0));
    Real n_win = fmax(dot(frame.n,dir_in),Real(0));
    // Real n_wout = fmax(dot(frame.n,dir_out),Real(0));
    Spectrum fm = schlick_fresnel(baseColor,h_wout);
    Real dm = GGX(frame,half,anistropic,roughness);
    Vector3 vlocal_in = to_local(frame,dir_in);
    Vector3 vlocal_out = to_local(frame,dir_out);
    Real gmw_in = smith_masking_gtr2(vlocal_in,roughness,anistropic);
    Real gmw_out = smith_masking_gtr2(vlocal_out,roughness,anistropic);
    Real gm = gmw_in*gmw_out;
    Spectrum metal = fm*dm*gm / (4*n_win);
    return metal;
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anistropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    
    // We use the reflectance to determine whether to choose specular sampling lobe or diffuse.
    // Real spec_prob = lS / (lS + lR);
    // Real diff_prob = 1 - spec_prob;
   
    // For the specular lobe, we use the ellipsoidal sampling from Heitz 2018
    // "Sampling the GGX Distribution of Visible Normals"
    // https://jcgt.org/published/0007/04/01/
    // this importance samples smith_masking(cos_theta_in) * GTR2(cos_theta_h, roughness) * cos_theta_out
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    if (n_dot_out <= 0 || n_dot_h <= 0) {
        return 0;
    }
    Real G = smith_masking_gtr2(to_local(frame, dir_in), roughness,anistropic);
    Real D = GGX(frame,half_vector,anistropic,roughness);
    // (4 * cos_theta_v) is the Jacobian of the reflectiokn
    Real spec_prob = (G * D) / (4 * n_dot_in);
    
    // For the diffuse lobe, we importance sample cos_theta_out
    // diff_prob *= n_dot_out / c_PI;
    // return spec_prob + diff_prob;
    return spec_prob;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
// Homework 1: implement this!
// Convert the incoming direction to local coordinates
    Vector3 local_dir_in = to_local(frame, dir_in);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real alpha = roughness * roughness;
    Vector3 local_micro_normal =
        sample_visible_normals(local_dir_in, alpha, rnd_param_uv);
    
    // Transform the micro normal to world space
    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, roughness /* roughness */
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
