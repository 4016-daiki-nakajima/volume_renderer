#include "../materials/microfacet.h"

Spectrum eval_op::operator()(const DisneySheen &bsdf) const {
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

    Real n_dot_out = dot(frame.n, dir_out);
    if (n_dot_out <= 0) return make_zero_spectrum();

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real h_dot_out      = dot(half_vector, dir_out);

    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool); 
    Real sheen_tint     = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real lum         = luminance(base_color);
    Spectrum C_tint  = lum > 0 ? base_color / lum : make_const_spectrum(1);
    Spectrum C_sheen = (1 - sheen_tint) + sheen_tint * C_tint;

    return C_sheen * pow(1 - fabs(h_dot_out), 5) * n_dot_out;
}

Real pdf_sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
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

    // cosine hemisphere importance sampling
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI; 
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // cosine hemisphere importance sampling
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0), Real(1) };
}

TextureSpectrum get_texture_op::operator()(const DisneySheen &bsdf) const {
    return bsdf.base_color;
}
