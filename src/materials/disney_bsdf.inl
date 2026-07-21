#include "../materials/microfacet.h"

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }

    // use eta if going into the surface, otherwise 1/eta
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    Spectrum base_color  = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_trans  = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool); 
    Real metallic        = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface      = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular        = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool); 
    Real roughness       = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint   = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic     = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool); 
    Real sheen           = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint      = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat       = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);  
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    
    // Clamp roughness to avoid numerical issues
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Real alpha_x;
    Real alpha_y;
    disney_aniso_alphas(roughness, anisotropic, alpha_x, alpha_y);

    // glass lobe
    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        half_vector = normalize(dir_in + dir_out * eta);
    }
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    Real h_dot_in  = dot(half_vector, dir_in); 
    Real h_dot_out = dot(half_vector, dir_out);
    
    Real F_g = fresnel_dielectric(h_dot_in, eta);
    Real D_g = GTR2_aniso(to_local(frame, half_vector), alpha_x, alpha_y);
    Real G_g = smith_masking_aniso(to_local(frame, dir_in), alpha_x, alpha_y) *
               smith_masking_aniso(to_local(frame, dir_out), alpha_x, alpha_y);

    Spectrum f_glass;
    if (reflect) {
        f_glass = base_color * (F_g * D_g * G_g) / (4 * fabs(dot(frame.n, dir_in))); 
    } else {
        Real eta_factor = dir == TransportDirection::TO_LIGHT ? (1 / (eta * eta)) : 1;
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Spectrum sqrt_base_color{sqrt(base_color.x), sqrt(base_color.y), sqrt(base_color.z)};
        f_glass = sqrt_base_color *
            (eta_factor * (1 - F_g) * D_g * G_g * eta * eta * fabs(h_dot_out * h_dot_in)) /
            (fabs(dot(frame.n, dir_in)) * sqrt_denom * sqrt_denom);
    }
    Real glass_weight = (1 - metallic) * specular_trans;
    Spectrum result = glass_weight * f_glass;

    // Diffuse, sheen, metal, clearcoat
    if (!reflect) return result;

    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    // if ray comes from inside the object, only the glass remains
    bool front_facing = dot(vertex.geometric_normal, dir_in) > 0;
    if (!front_facing || n_dot_out <= 0) return result;

    // ---- Diffuse ----
    Real Fd90   = Real(0.5) + 2 * roughness * h_dot_out * h_dot_out;
    Real FD_in  = 1 + (Fd90 - 1) * pow(1 - fabs(n_dot_in), 5);
    Real FD_out = 1 + (Fd90 - 1) * pow(1 - fabs(n_dot_out), 5);
    Spectrum base_diffuse = (base_color / c_PI) * FD_in * FD_out * n_dot_out;

    Real Fss90   = roughness * h_dot_out * h_dot_out;
    Real Fss_in  = 1 + (Fss90 - 1) * pow(1 - fabs(n_dot_in), 5);
    Real Fss_out = 1 + (Fss90 - 1) * pow(1 - fabs(n_dot_out), 5);
    Spectrum subsurface_term = (Real(1.25) * base_color / c_PI) *
        (Fss_in * Fss_out * (1 / (fabs(n_dot_in) + fabs(n_dot_out)) - Real(0.5)) + Real(0.5)) *
        n_dot_out;

    Spectrum f_diffuse = (1 - subsurface) * base_diffuse + subsurface * subsurface_term;

    // ---- Sheen ----
    Real lum = luminance(base_color);
    Spectrum C_tint = lum > 0 ? base_color / lum : make_const_spectrum(1);
    Spectrum C_sheen = (1 - sheen_tint) + sheen_tint * C_tint;
    Spectrum f_sheen = C_sheen * pow(1 - fabs(h_dot_out), 5) * n_dot_out;

    // ---- Metal with dielectric term ----
    Real R0_eta  = (eta - 1) * (eta - 1) / ((eta + 1) * (eta + 1));
    Spectrum Ks  = (1 - specular_tint) + specular_tint * C_tint;
    Spectrum C0  = specular * R0_eta * (1 - metallic) * Ks + metallic * base_color;
    Spectrum F_m = schlick_fresnel(C0, h_dot_out);
    Spectrum f_metal = F_m * D_g * G_g / (4 * n_dot_in);

    // ---- Clearcoat ----
    Real alpha_g = clearcoat_alpha_g(clearcoat_gloss);
    Real n_dot_h = dot(frame.n, half_vector);
    Real Fc = schlick_fresnel(Real(0.04), h_dot_out);
    Real Dc = GTR1(n_dot_h, alpha_g);
    Real Gc = smith_masking_clearcoat(to_local(frame, dir_in)) *
              smith_masking_clearcoat(to_local(frame, dir_out));
    Real f_clearcoat = (Fc * Dc * Gc) / (4 * n_dot_in);

    result +=
        (1 - specular_trans) * (1 - metallic) * f_diffuse +
        (1 - metallic) * sheen * f_sheen +
        (1 - specular_trans * (1 - metallic)) * f_metal +
        Real(0.25) * clearcoat * make_const_spectrum(f_clearcoat);

    return result;
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }

    // use eta if going into the surface, otherwise 1/eta
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    Real specular_trans  = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic        = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool); 
    Real roughness       = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool); 
    Real anisotropic     = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat       = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);  
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Real alpha_x, alpha_y;
    disney_aniso_alphas(roughness, anisotropic, alpha_x, alpha_y);

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        half_vector = normalize(dir_in + dir_out * eta);
    }
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }
    Real h_dot_in = dot(half_vector, dir_in);
    Real h_dot_out = dot(half_vector, dir_out);

    Real F_g  = fresnel_dielectric(h_dot_in, eta); 
    Real D_g  = GTR2_aniso(to_local(frame, half_vector), alpha_x, alpha_y);
    Real G_in = smith_masking_aniso(to_local(frame, dir_in), alpha_x, alpha_y);

    // PDF of the glass lobe's reflect/refract sampling
    Real pdf_glass;
    if (reflect) {
        pdf_glass = (F_g * D_g * G_in) / (4 * fabs(dot(frame.n, dir_in)));
    } else {
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom); 
        pdf_glass = (1 - F_g) * D_g * G_in * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
    }

    // selection weights for the four sampling strategies
    bool front_facing     = dot(vertex.geometric_normal, dir_in) > 0;
    Real diffuse_weight   = (1 - metallic) * (1 - specular_trans);
    Real metal_weight     = (1 - specular_trans * (1 - metallic));
    Real glass_weight     = (1 - metallic) * specular_trans;
    Real clearcoat_weight = Real(0.25) * clearcoat;
    if (!front_facing) {
        // if the ray comes from inside the object, only the glass lobe remains
        diffuse_weight = metal_weight = clearcoat_weight = 0;
        glass_weight = 1;
    }
    Real total_weight = diffuse_weight + metal_weight + glass_weight + clearcoat_weight;
    
    if (total_weight <= 0) return 0;
     
    Real p_diffuse   = diffuse_weight / total_weight; 
    Real p_metal     = metal_weight / total_weight;    
    Real p_glass     = glass_weight / total_weight;
    Real p_clearcoat = clearcoat_weight / total_weight;  

    Real pdf = p_glass * pdf_glass;
    if (reflect) {
        Real n_dot_in = dot(frame.n, dir_in);
        Real n_dot_out = dot(frame.n, dir_out);
        if (n_dot_out > 0) {
            pdf += p_diffuse * (n_dot_out / c_PI);
            pdf += p_metal * (G_in * D_g) / (4 * n_dot_in);

            Real n_dot_h = dot(frame.n, half_vector);
            Real alpha_g = clearcoat_alpha_g(clearcoat_gloss);
            Real Dc = GTR1(n_dot_h, alpha_g);
            pdf += p_clearcoat * (Dc * n_dot_h) / (4 * fabs(h_dot_out));
        }
    }
    return pdf;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // If going into the surface, then we use normal eta (internal/external), otherwise external/internal
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }

    Real specular_trans  = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);  
    Real metallic        = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool); 
    Real roughness       = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic     = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);  
    Real clearcoat       = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);  
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    // Selection weights for the four sampling strategies (unnormalized).
    bool front_facing     = dot(vertex.geometric_normal, dir_in) > 0;
    Real diffuse_weight   = (1 - metallic) * (1 - specular_trans);
    Real metal_weight     = (1 - specular_trans * (1 - metallic));
    Real glass_weight     = (1 - metallic) * specular_trans;
    Real clearcoat_weight = Real(0.25) * clearcoat;
    if (!front_facing) {
        // if the ray comes from inside the object, only the glass lobe remains
        diffuse_weight = metal_weight = clearcoat_weight = 0;
        glass_weight = 1;
    }
    Real total_weight = diffuse_weight + metal_weight + glass_weight + clearcoat_weight;
    if (total_weight <= 0) {
        return {};
    }
    Real p_diffuse = diffuse_weight / total_weight;
    Real p_metal   = metal_weight / total_weight;
    Real p_glass   = glass_weight / total_weight;

    Real alpha_x, alpha_y;
    disney_aniso_alphas(roughness, anisotropic, alpha_x, alpha_y);

    Real w = rnd_param_w;
    if (w < p_diffuse) {
        // ---- Diffuse lobe ----
        return BSDFSampleRecord{
            to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
            Real(0), roughness};
    } else if (w < p_diffuse + p_metal) {
        // ---- Metal lobe ----
        Vector3 local_dir_in = to_local(frame, dir_in);
        Vector3 local_micro_normal =
            sample_visible_normals_aniso(local_dir_in, alpha_x, alpha_y, rnd_param_uv);
        Vector3 half_vector = to_world(frame, local_micro_normal);
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        return BSDFSampleRecord{reflected, Real(0), roughness};
    } else if (w < p_diffuse + p_metal + p_glass) {
        // ---- Glass lobe ----
        Real w_glass = (w - (p_diffuse + p_metal)) / p_glass;

        Vector3 local_dir_in = to_local(frame, dir_in);
        Vector3 local_micro_normal =
            sample_visible_normals_aniso(local_dir_in, alpha_x, alpha_y, rnd_param_uv);
        Vector3 half_vector = to_world(frame, local_micro_normal);
        if (dot(half_vector, frame.n) < 0) half_vector = -half_vector;

        Real h_dot_in = dot(half_vector, dir_in);
        Real F = fresnel_dielectric(h_dot_in, eta);

        if (w_glass <= F) {
            // Reflection
            Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
            return BSDFSampleRecord{reflected, Real(0), roughness};
        } else {
            // Refraction
            Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
            if (h_dot_out_sq <= 0)  return {};
            if (h_dot_in < 0) half_vector = -half_vector;

            Real h_dot_out    = sqrt(h_dot_out_sq);
            Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
            return BSDFSampleRecord{refracted, eta, roughness};
        }
    } else {
        // ---- Clearcoat ----
        Real alpha_g = clearcoat_alpha_g(clearcoat_gloss);
        Real alpha2 = alpha_g * alpha_g;
        Real u0 = rnd_param_uv.x, u1 = rnd_param_uv.y;
        Real cos_h_elevation = sqrt(std::max(Real(0), (1 - pow(alpha2, 1 - u0)) / (1 - alpha2)));
        Real sin_h_elevation = sqrt(std::max(Real(0), 1 - cos_h_elevation * cos_h_elevation));
        Real h_azimuth = 2 * c_PI * u1;
        Vector3 local_micro_normal{
            sin_h_elevation * cos(h_azimuth),
            sin_h_elevation * sin(h_azimuth),
            cos_h_elevation};

        Vector3 half_vector = to_world(frame, local_micro_normal);
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        return BSDFSampleRecord{reflected, Real(0), alpha_g};
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
