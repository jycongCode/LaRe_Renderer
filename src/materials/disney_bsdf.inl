#include "../microfacet.h"
#include "../table_dist.h"

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    
    Spectrum fglass;
    {
        Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
        Real anistropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
        
        Spectrum Ks = eval(
            bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
        Spectrum Kt = Spectrum(sqrt(Ks.x),sqrt(Ks.y),sqrt(Ks.z));
        Real roughness = eval(
            bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);

        Vector3 half_vector;
        if (reflect) {
            half_vector = normalize(dir_in + dir_out);
        } else {
            // "Generalized half-vector" from Walter et al.
            // See "Microfacet Models for Refraction through Rough Surfaces"
            half_vector = normalize(dir_in + dir_out * eta);
        }

        // Flip half-vector if it's below surface
        if (dot(half_vector, frame.n) < 0) {
            half_vector = -half_vector;
        }

        // Clamp roughness to avoid numerical issues.
        roughness = std::clamp(roughness, Real(0.01), Real(1));

        // Compute F / D / G
        Real h_dot_in = dot(half_vector, dir_in);
        Real F = fresnel_dielectric(h_dot_in, eta);
        Real D = GGX(frame,half_vector,anistropic,roughness);
        Vector3 vlocal_in = to_local(frame,dir_in);
        Vector3 vlocal_out = to_local(frame,dir_out);
        Real gmw_in = smith_masking_gtr2(vlocal_in,roughness,anistropic);
        Real gmw_out = smith_masking_gtr2(vlocal_out,roughness,anistropic);
        Real G = gmw_in*gmw_out;
        if (reflect) {
            fglass = Ks * (F * D * G) / (4 * fabs(dot(frame.n, dir_in)));
        } else {
            Real eta_factor = dir == TransportDirection::TO_LIGHT ? (1 / (eta * eta)) : 1;
            Real h_dot_out = dot(half_vector, dir_out);
            Real sqrt_denom = h_dot_in + eta * h_dot_out;
           
            fglass = Kt * (eta_factor * (1 - F) * D * G * eta * eta * fabs(h_dot_out * h_dot_in)) / 
                (fabs(dot(frame.n, dir_in)) * sqrt_denom * sqrt_denom);
        }
    }

    // Homework 1: implement this!
    Spectrum fdiffuse;
    {
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
        
        fdiffuse = (1.0-subsurface)*diffuse + subsurface*subsf;
    }

    Spectrum fmetal;
    {
        // sample texture
        Real anistropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
        Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
        
        // vector precompute
        Vector3 half = normalize(dir_in + dir_out);
        Real h_wout = fmax(dot(half, dir_out), Real(0));
        Real n_win = fmax(dot(frame.n,dir_in),Real(0));

        // fm with specular tint
        Spectrum c_tint = base_color * ( luminance(base_color)>0?(Real(1)/luminance(base_color)):Real(1));
        Spectrum mKs = (Real(1)-specular_tint) + specular_tint * c_tint;
        Spectrum mC0 = specular * R0(bsdf.eta)*(Real(1)-metallic)*mKs + metallic * base_color;
        Spectrum fm = schlick_fresnel(mC0,h_wout);
        
        // original metal dm gm
        Real dm = GGX(frame,half,anistropic,roughness);
        Vector3 vlocal_in = to_local(frame,dir_in);
        Vector3 vlocal_out = to_local(frame,dir_out);
        Real gmw_in = smith_masking_gtr2(vlocal_in,roughness,anistropic);
        Real gmw_out = smith_masking_gtr2(vlocal_out,roughness,anistropic);
        Real gm = gmw_in*gmw_out;
        fmetal = fm*dm*gm / (4*n_win);
    }

    Spectrum fsheen;
    {
        Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
        Spectrum base_color = eval(bsdf.base_color,vertex.uv,vertex.uv_screen_size,texture_pool);
        Vector3 half_vector = normalize(dir_in + dir_out);

        Spectrum c_tint = base_color * ((Real(1)/luminance(base_color)) ? luminance(base_color)>0:Real(1));
        Spectrum c_sheen = (Real(1)-sheen_tint)+sheen_tint*c_tint;

        Real n_wout = fmax(dot(frame.n,dir_out),Real(0));
        Real h_wout = fmax(dot(half_vector,dir_out),Real(0));

        fsheen = c_sheen * pow(Real(1)-h_wout,5)*n_wout;
    }

    Spectrum fclearcoat;
    {
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
        Spectrum(clearcoat,clearcoat,clearcoat);
    }

    // texture sampling
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);

    Spectrum fdisney;
    if(dot(dir_in,vertex.geometric_normal)>Real(0)){
        fdisney = (Real(1)-specular_transmission)*(Real(1)-metallic)*fdiffuse+\
                        (Real(1)-metallic)*sheen*fsheen+\
                        (Real(1)-specular_transmission*(Real(1)-metallic))*fmetal+\
                        0.25*clearcoat*fclearcoat+\
                        (Real(1)-metallic)*specular_transmission*fglass;
    }else{
        fdisney = fglass;
    }
    // fdisney = fglass;
    return fdisney;
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    
    Real fglass;
    {
        // Homework 1: implement this!
        // If we are going into the surface, then we use normal eta
        // (internal/external), otherwise we use external/internal.
        Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
        assert(eta > 0);

        Vector3 half_vector;
        if (reflect) {
            half_vector = normalize(dir_in + dir_out);
        } else {
            // "Generalized half-vector" from Walter et al.
            // See "Microfacet Models for Refraction through Rough Surfaces"
            half_vector = normalize(dir_in + dir_out * eta);
        }

        // Flip half-vector if it's below surface
        if (dot(half_vector, frame.n) < 0) {
            half_vector = -half_vector;
        }

        Real roughness = eval(
            bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
        // Clamp roughness to avoid numerical issues.
        roughness = std::clamp(roughness, Real(0.01), Real(1));

        // We sample the visible normals, also we use F to determine
        // whether to sample reflection or refraction
        // so PDF ~ F * D * G_in for reflection, PDF ~ (1 - F) * D * G_in for refraction.
        Real h_dot_in = dot(half_vector, dir_in);
        Real F = fresnel_dielectric(h_dot_in, eta);
        Real anistropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real D = GGX(frame,half_vector,anistropic,roughness);
        Vector3 vlocal_in = to_local(frame,dir_in);
        Vector3 vlocal_out = to_local(frame,dir_out);
        Real G_in = smith_masking_gtr2(vlocal_in,roughness,anistropic);
        
        if (reflect) {
            fglass = (F * D * G_in) / (4 * fabs(dot(frame.n, dir_in)));
        } else {
            Real h_dot_out = dot(half_vector, dir_out);
            Real sqrt_denom = h_dot_in + eta * h_dot_out;
            Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
            fglass = (1 - F) * D * G_in * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
        }
    }
    
    // diffuse pdf
    Real fdiffuse = fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
   
    // metal specular pdf
    Real fmetal = 0.0;
    {
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
            fmetal = 0.0;
        }
        Real G = smith_masking_gtr2(to_local(frame, dir_in), roughness,anistropic);
        Real D = GGX(frame,half_vector,anistropic,roughness);
        // (4 * cos_theta_v) is the Jacobian of the reflectiokn
        fmetal = (G * D) / (4 * n_dot_in);
    }
    
    // fclearcoat
    Real fclearcoat = 0.0;
    {
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
        fclearcoat = dc*n_h / (Real(4)*h_wout);
    }

    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);

    // Homework 1: implement this!
    Real diffuse_weight = (Real(1)-specular_transmission)*(Real(1)-metallic);
    Real metal_weight =  (Real(1)-specular_transmission*(Real(1)-metallic));
    Real glass_weight = (Real(1)-metallic)*specular_transmission;
    Real clearcoat_weight = 0.25*clearcoat;

    std::vector<Real> dist = {glass_weight,diffuse_weight,metal_weight,clearcoat_weight};
    TableDist1D table1d = make_table_dist_1d(dist);
    Real fdisney;

    if(dot(dir_in,vertex.geometric_normal) > Real(0)){
        fdisney = table1d.pmf[1] * fdiffuse+\
                    table1d.pmf[2] * fmetal+\
                     table1d.pmf[0] * fglass+\
                      table1d.pmf[3] * fclearcoat;
    }else{
        fdisney = fglass;
    }
    // fdisney = fglass;
    return fdisney;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    
    Real diffuse_weight = (Real(1)-specular_transmission)*(Real(1)-metallic);
    Real metal_weight =  (Real(1)-specular_transmission*(Real(1)-metallic));
    Real glass_weight = (Real(1)-metallic)*specular_transmission;
    Real clearcoat_weight = 0.25*clearcoat;

    std::vector<Real> dist = {glass_weight,diffuse_weight,metal_weight,clearcoat_weight};
    TableDist1D table1d = make_table_dist_1d(dist);
    int idx = sample(table1d,rnd_param_w);
    // if ray comes from inside, apply glass lobe
    if(dot(dir_in,vertex.geometric_normal)<=0)idx = 0;
    // idx=0;
    if(idx == 1){
        return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, Real(1) /* roughness */};
    }else if(idx == 2){
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
    }else if(idx == 0){
        Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
        // Homework 1: implement this!
        Real roughness = eval(
            bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
        // Clamp roughness to avoid numerical issues.
        roughness = std::clamp(roughness, Real(0.01), Real(1));
        // Sample a micro normal and transform it to world space -- this is our half-vector.
        Real alpha = roughness * roughness;
        Vector3 local_dir_in = to_local(frame, dir_in);
        Vector3 local_micro_normal =
            sample_visible_normals(local_dir_in, alpha, rnd_param_uv);

        Vector3 half_vector = to_world(frame, local_micro_normal);
        // Flip half-vector if it's below surface
        if (dot(half_vector, frame.n) < 0) {
            half_vector = -half_vector;
        }

        // Now we need to decide whether to reflect or refract.
        // We do this using the Fresnel term.
        Real h_dot_in = dot(half_vector, dir_in);
        Real F = fresnel_dielectric(h_dot_in, eta);

        Real remap_rnd_param;
        if(rnd_param_w < table1d.pmf[0]){
            remap_rnd_param = rnd_param_w / table1d.pmf[0];
        }else{
            remap_rnd_param = (rnd_param_w-table1d.pmf[0]) / (Real(1)-table1d.pmf[0]);
        }
        if (remap_rnd_param <= F) {
            // Reflection
            Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
            // set eta to 0 since we are not transmitting
            return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
        } else {
            // Refraction
            // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
            // (note that our eta is eta2 / eta1, and l = -dir_in)
            Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
            if (h_dot_out_sq <= 0) {
                // Total internal reflection
                // This shouldn't really happen, as F will be 1 in this case.
                return {};
            }
            // flip half_vector if needed
            if (h_dot_in < 0) {
                half_vector = -half_vector;
            }
            Real h_dot_out= sqrt(h_dot_out_sq);
            Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
            return BSDFSampleRecord{refracted, eta, roughness};
        }
    }else if(idx == 3){
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
    return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
