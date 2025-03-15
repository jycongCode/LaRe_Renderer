#include "../microfacet.h"
Real R0(Real mu=1.5){
    return (mu-1.0)*(mu-1.0) / ((mu+1.0)*(mu+1.0));
}

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3 half = normalize(dir_in+dir_out);
    Real h_wout = fmax(dot(half,dir_out),Real(0));
    Real n_win = fmax(dot(frame.n,dir_in),Real(0));
    Vector3 h_local = to_local(frame,half);
    
    Real fc = schlick_fresnel(R0(1.5),h_wout);
    Real alphag = (Real(1)-clearcoat_gloss)*0.1+clearcoat_gloss*0.001;
    Real ag2 = alphag * alphag;
    Real dc = (ag2-1)/(c_PI*log(ag2)*(Real(1)+(ag2-1)*h_local.z*h_local.z));
    Real gw_in = smith_masking_gtr2(to_local(frame,dir_in),0.25);
    Real gw_out = smith_masking_gtr2(to_local(frame,dir_out),0.25);
    Real gc = gw_in*gw_out;
    Real clearcoat = fc*dc*gc/(Real(4)*n_win);
    return Spectrum(clearcoat,clearcoat,clearcoat);
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3 half = normalize(dir_in+dir_out);
    Real h_wout = fmax(dot(half,dir_out),Real(0));
    Real n_win = fmax(dot(frame.n,dir_in),Real(0));
    Real n_h = fmax(dot(half,frame.n),Real(0));
    Vector3 h_local = to_local(frame,half);
    
    Real fc = schlick_fresnel(R0(1.5),h_wout);
    Real alphag = (Real(1)-clearcoat_gloss)*0.1+clearcoat_gloss*0.001;
    Real ag2 = alphag * alphag;
    Real dc = (ag2-1)/(c_PI*log(ag2)*(Real(1)+(ag2-1)*h_local.z*h_local.z));
    return dc*n_h / (Real(4)*h_wout);
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha = (Real(1)-clearcoat_gloss)*0.1+clearcoat_gloss*0.001;
    Real a2 = alpha*alpha;
    Real cos_he = sqrt((Real(1)-pow(a2,Real(1-rnd_param_uv.x)))/(Real(1)-a2));
    Real ha = Real(2)*c_PI*rnd_param_uv.y;
    Real sin_he = sqrt(fmax(1-cos_he*cos_he,0.0));
    Vector3 h_local(sin_he*cos(ha),sin_he*sin(ha),cos_he);
    Vector3 half = to_world(frame,h_local);
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half) * half);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, Real(1) /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
