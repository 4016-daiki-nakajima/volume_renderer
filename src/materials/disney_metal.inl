#include "../materials/microfacet.h"

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

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_in  = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out); 
    Real n_dot_h   = dot(frame.n, half_vector);
    if (n_dot_out <= 0 || n_dot_h <= 0) {
        return make_zero_spectrum();
    }

    Spectrum base_color = eval(bsdf.base_color,  vertex.uv, vertex.uv_screen_size, texture_pool); 
    Real roughness      = eval(bsdf.roughness,   vertex.uv, vertex.uv_screen_size, texture_pool); 
    Real anisotropic    = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool); 
    
    // clamp roughness to avoid numerical issues 
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Real alpha_x, alpha_y;
    disney_aniso_alphas(roughness, anisotropic, alpha_x, alpha_y);

    Real h_dot_out = dot(half_vector, dir_out);
    Spectrum F_m = schlick_fresnel(base_color, h_dot_out);
    Real D_m = GTR2_aniso(to_local(frame, half_vector), alpha_x, alpha_y);
    Real G_m = smith_masking_aniso(to_local(frame, dir_in), alpha_x, alpha_y) *
               smith_masking_aniso(to_local(frame, dir_out), alpha_x, alpha_y);

    return F_m * D_m * G_m / (4 * n_dot_in);
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

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_in = dot(frame.n, dir_in);  
    Real n_dot_h  = dot(frame.n, half_vector);
    if (dot(frame.n, dir_out) <= 0 || n_dot_h <= 0) {
        return 0;
    }
 
    Real roughness   = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Real alpha_x, alpha_y; 
    disney_aniso_alphas(roughness, anisotropic, alpha_x, alpha_y);
 
    Real G_in = smith_masking_aniso(to_local(frame, dir_in), alpha_x, alpha_y);
    Real D = GTR2_aniso(to_local(frame, half_vector), alpha_x, alpha_y);
    return (G_in * D) / (4 * n_dot_in);
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

    Real roughness   = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv,  vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Real alpha_x, alpha_y;
    disney_aniso_alphas(roughness, anisotropic, alpha_x, alpha_y);

    Vector3 local_dir_in = to_local(frame, dir_in);  
    Vector3 local_micro_normal =
        sample_visible_normals_aniso(local_dir_in, alpha_x, alpha_y, rnd_param_uv);

    // Transform micro normal to world space  
    Vector3 half_vector = to_world(frame, local_micro_normal); 
    // Reflect over the world space half vector 
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        Real(0), roughness };
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
