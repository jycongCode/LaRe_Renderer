inline Real quickpow(Real x,int y){
    Real ans = 1.0;
    while(y){
        if(y&1)ans*=x;
        x*=x;
        y >>= 1;
    }
    return ans;
}

inline Real FD(Vector3 dir,Vector3 normal,Real fd90){
    Real n_w = fmax(dot(normal,dir),Real(0));
    return (Real(1)+(fd90-Real(1))*quickpow(1-n_w,5));
}

inline Real FSS(Vector3 dir,Vector3 normal,Real fss90){
    Real n_w = fmax(dot(normal,dir),Real(0));
    return (Real(1)+(fss90-Real(1))*quickpow(1-n_w,5));
}

Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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

    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum baseColor = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real h_wout = fmax(dot(half_vector, dir_out), Real(0));
    Real n_win = fmax(dot(frame.n,dir_in),Real(0));
    Real n_wout = fmax(dot(frame.n,dir_out),Real(0));
    
    Real baseFD90 = Real(0.5) + Real(2.0)*roughness*h_wout*h_wout;
    Real baseFD_win = FD(dir_in,frame.n,baseFD90);
    Real baseFD_wout = FD(dir_out,frame.n,baseFD90);
    Spectrum diffuse = baseColor * Real(1.0/c_PI)*baseFD_win*baseFD_wout*n_wout;
    
    Real subFss90 = roughness*h_wout*h_wout;
    Real subFss_win = FSS(dir_in,frame.n,subFss90);
    Real subFss_wout = FSS(dir_out,frame.n,subFss90);
    Spectrum subsf = baseColor*Real(1.25/c_PI)*(subFss_win*subFss_wout*(Real(1.0/(n_win+n_wout))-0.5)+0.5)*n_wout;
    
    return (1.0-subsurface)*diffuse + subsurface*subsf;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
     return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, Real(1) /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
