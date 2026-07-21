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

    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    if (n_dot_out <= 0) {
        return make_zero_spectrum();
    }

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real h_dot_out = dot(half_vector, dir_out);

    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);

    // Base diffuse lobe
    Real Fd90   = Real(0.5) + 2 * roughness * h_dot_out * h_dot_out;
    Real FD_in  = 1 + (Fd90 - 1) * pow(1 - fabs(n_dot_in), 5);
    Real FD_out = 1 + (Fd90 - 1) * pow(1 - fabs(n_dot_out), 5);
    Spectrum base_diffuse = (base_color / c_PI) * FD_in * FD_out * n_dot_out;

    // Subsurface lobe
    Real Fss90   = roughness * h_dot_out * h_dot_out;
    Real Fss_in  = 1 + (Fss90 - 1) * pow(1 - fabs(n_dot_in), 5);
    Real Fss_out = 1 + (Fss90 - 1) * pow(1 - fabs(n_dot_out), 5); 
    Spectrum subsurface_term = (Real(1.25) * base_color / c_PI) *
        (Fss_in * Fss_out * (1 / (fabs(n_dot_in) + fabs(n_dot_out)) - Real(0.5)) + Real(0.5)) *
        n_dot_out;

    return (1 - subsurface) * base_diffuse + subsurface * subsurface_term;
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

    // For Disney diffuse, we importance sample the cosine hemisphere domain.
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

    // cosine hemisphere importance sampling
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0), roughness };
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
