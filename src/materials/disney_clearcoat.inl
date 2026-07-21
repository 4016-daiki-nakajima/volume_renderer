#include "../materials/microfacet.h"

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

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    if (n_dot_out <= 0 || n_dot_h <= 0) {
        return make_zero_spectrum();
    }

    Real gloss   = eval( bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_g = clearcoat_alpha_g(gloss);

    // R0(eta = 1.5) = ((1.5 - 1) / (1.5 + 1))^2 = 0.04
    Real R0 = Real(0.04);
    Real h_dot_out = dot(half_vector, dir_out);
    Real Fc = schlick_fresnel(R0, h_dot_out);
    Real Dc = GTR1(n_dot_h, alpha_g);
    Real Gc = smith_masking_clearcoat(to_local(frame, dir_in)) *
              smith_masking_clearcoat(to_local(frame, dir_out));

    Real f = (Fc * Dc * Gc) / (4 * n_dot_in);
    return make_const_spectrum(f);
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

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    if (dot(frame.n, dir_out) <= 0 || n_dot_h <= 0)
        return 0;

    Real gloss   = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_g = clearcoat_alpha_g(gloss);

    Real Dc = GTR1(n_dot_h, alpha_g);
    Real h_dot_out = dot(half_vector, dir_out);

    return (Dc * n_dot_h) / (4 * fabs(h_dot_out));
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

    Real gloss   = eval( bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_g = clearcoat_alpha_g(gloss);
    Real alpha2  = alpha_g * alpha_g;

    Real u0 = rnd_param_uv.x;
    Real u1 = rnd_param_uv.y;
    Real cos_h_elev = sqrt(std::max(Real(0), (1 - pow(alpha2, 1 - u0)) / (1 - alpha2)));
    Real sin_h_elev = sqrt(std::max(Real(0), 1 - cos_h_elev * cos_h_elev));
    Real h_az = 2 * c_PI * u1;
    Vector3 local_micro_normal{
        sin_h_elev * cos(h_az),
        sin_h_elev * sin(h_az),
        cos_h_elev};

    // Transform micro normal to world space
    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Reflect over  world space half vector
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{ reflected, Real(0), alpha_g };
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
