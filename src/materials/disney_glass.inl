#include "../materials/microfacet.h"

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // If going into the surface, then we use normal eta  (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Real alpha_x, alpha_y;
    disney_aniso_alphas(roughness, anisotropic, alpha_x, alpha_y);

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        half_vector = normalize(dir_in + dir_out * eta);
    }
    // Flip half-vector if below surface
    if (dot(half_vector, frame.n) < 0) half_vector = -half_vector; 

    Real h_dot_in = dot(half_vector, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);
    Real D = GTR2_aniso(to_local(frame, half_vector), alpha_x, alpha_y);
    Real G = smith_masking_aniso(to_local(frame, dir_in), alpha_x, alpha_y) *
             smith_masking_aniso(to_local(frame, dir_out), alpha_x, alpha_y);

    if (reflect) {
        return base_color * (F * D * G) / (4 * fabs(dot(frame.n, dir_in)));
    } else {
        Real eta_factor = dir == TransportDirection::TO_LIGHT ? (1 / (eta * eta)) : 1;
        Real h_dot_out  = dot(half_vector, dir_out);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Spectrum sqrt_base_color{sqrt(base_color.x), sqrt(base_color.y), sqrt(base_color.z)};
        
        return sqrt_base_color * (eta_factor * (1 - F) * D * G * eta * eta * fabs(h_dot_out * h_dot_in)) /
            (fabs(dot(frame.n, dir_in)) * sqrt_denom * sqrt_denom);
    }
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // If going into the surface, then we use normal eta (internal/external), otherwise external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        half_vector = normalize(dir_in + dir_out * eta);
    }
    // Flip half-vector if below surface
    if (dot(half_vector, frame.n) < 0) half_vector = -half_vector; 

    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Real alpha_x, alpha_y;
    disney_aniso_alphas(roughness, anisotropic, alpha_x, alpha_y);

    // compute the PDF of the half-vector.
    Real h_dot_in = dot(half_vector, dir_in);
    Real F        = fresnel_dielectric(h_dot_in, eta);
    Real D         = GTR2_aniso(to_local(frame, half_vector), alpha_x, alpha_y);
    Real G_in      = smith_masking_aniso(to_local(frame, dir_in), alpha_x, alpha_y);
    if (reflect) {
        return (F * D * G_in) / (4 * fabs(dot(frame.n, dir_in)));
    } else {
        Real h_dot_out  = dot(half_vector, dir_out);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real dh_dout    = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
        return (1 - F) * D * G_in * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
    }
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    // If going into the surface, then we use normal eta (internal/external), otherwise external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }

    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues
    roughness = std::clamp(roughness, Real(0.01), Real(1)); 

    Real alpha_x, alpha_y;
    disney_aniso_alphas(roughness, anisotropic, alpha_x, alpha_y);

    // Sample micro normal and transform it to world space
    Vector3 local_dir_in       = to_local(frame, dir_in);
    Vector3 local_micro_normal = sample_visible_normals_aniso(local_dir_in, alpha_x, alpha_y, rnd_param_uv);

    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Flip half-vector if below surface
    if (dot(half_vector, frame.n) < 0) half_vector = -half_vector; 

    // reflect or refract based on Fresnel
    Real h_dot_in = dot(half_vector, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);

    if (rnd_param_w <= F) {
        // Reflection
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        // Set eta to 0 since we are not transmitting.
        return BSDFSampleRecord{reflected, Real(0), roughness};
    } else {
        // Refraction
        Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
        if (h_dot_out_sq <= 0) {
            // total internal reflection (shouldn't happen)
            return {};
        }
        // Flip half_vector if needed.
        if (h_dot_in < 0) half_vector = -half_vector;
        Real h_dot_out    = sqrt(h_dot_out_sq);
        Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
        
        return BSDFSampleRecord{refracted, eta, roughness};
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
