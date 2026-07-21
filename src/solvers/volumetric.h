#pragma once


int update_medium(PathVertex intersection, Ray &ray, int medium) {
    int new_medium = medium;

    if (intersection.interior_medium_id != intersection.exterior_medium_id){
        // At medium transition. Update medium.
        if (dot(ray.dir, intersection.geometric_normal) > 0) {
            new_medium = intersection.exterior_medium_id;
        } else {
            new_medium = intersection.interior_medium_id;
        }
    }
    return new_medium;
}

Spectrum next_event_estimation(const Scene& scene, Vector3 p, Vector3 dir_view, int current_medium, 
                                 const Material *mat, const PathVertex *pathvertex, int bounces, pcg32_state &rng ){

    Real shape_rnd  = next_pcg32_real<Real>(rng); 
    Real light_rand = next_pcg32_real<Real>(rng);
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};

    int sampled_light  = sample_light(scene, light_rand);
    const Light &light = scene.lights[sampled_light];

    auto light_sample = sample_point_on_light(light, p, light_uv, shape_rnd, scene);
    Vector3 p_prime = light_sample.position;

    Vector3 dir_light = normalize(p_prime - p);
    Real dist_light   = distance(p, p_prime);
    Vector3 p_initial = p;

    // Compute transmittance to light. Skip through index-matching shapes.
    Spectrum T_light     = make_const_spectrum(1);
    int shadow_medium    = current_medium;
    int shadow_bounces   = 0;
    Spectrum p_trans_dir = make_const_spectrum(1); 
    Spectrum p_trans_nee = make_const_spectrum(1);



    while (true) {
        Real eps = get_shadow_epsilon(scene);
        Ray shadow_ray{p, dir_light, eps, (1 - eps) * dist_light};

        auto isect = intersect(scene, shadow_ray);
        Real next_t = distance(p, p_prime);

        if (isect) next_t = distance(p, isect->position); 
        
        // Account for the transmittance to next_t
        if (shadow_medium != -1) {
            const Medium &medium = scene.media[shadow_medium];
            Spectrum majorant = get_majorant(medium, shadow_ray);

            Real u = next_pcg32_real<Real>(rng);
            int channel = min(int(u * 3), 2);
            int iteration = 0;
            Real accum_t = 0;

            while (true) {
                if (majorant[channel] <= 0) break;
                if (iteration >= scene.options.max_null_collisions) break;

                Real t  = -log(1 - next_pcg32_real<Real>(rng)) / majorant[channel];
                Real dt = next_t - accum_t;
                accum_t = min(accum_t + t, next_t);
                if (t < dt) {
                    // didn't hit a surface, null-scattering event
                    p += + t * dir_light;
                    Spectrum sigma_t = get_sigma_a(medium, p) + get_sigma_s(medium, p);

                    T_light *= exp(-majorant * t) * (majorant - sigma_t) / max(majorant);
                    p_trans_nee *= exp(-majorant * t) * majorant / max(majorant);
                    Spectrum real_prob = sigma_t / majorant;
                    p_trans_dir *= exp(-majorant * t) * majorant * (1 - real_prob) / max(majorant);
                    if (max(T_light) <= 0) break;
                } else { // hit the surface
                    T_light     *= exp(-majorant * dt);
                    p_trans_nee *= exp(-majorant * dt);
                    p_trans_dir *= exp(-majorant * dt);
                    break;
                }
                iteration++;
            }
        }
        if (!isect) {
        // Nothing is blocking, we’re done
            break;
        } else {
            // Something is blocking: is it an opaque surface?
            if (isect->material_id >= 0) {
                // we're blocked
                return make_zero_spectrum();
            }
            // otherwise, it’s an index-matching surface and
            // we want to pass through -- this introduces
            // one extra connection vertex
            shadow_bounces++;
            if (scene.options.max_depth != -1 && bounces + shadow_bounces + 1 >= scene.options.max_depth) {
                // Reach the max no. of vertices
                return make_zero_spectrum();
            }
            shadow_medium = update_medium(*isect, shadow_ray, shadow_medium);
            p = isect->position;
        }
    }


    if (max(T_light) > 0) {

        // clamp negative values
        for (int i = 0; i < 3; i++) {
            if (T_light[i] <= 0) T_light[i] = 0;
        }

        if (mat) {
            if (!pathvertex) return make_zero_spectrum();

            Real G = max(-dot(dir_light, light_sample.normal), 0.0) / (dist_light * dist_light);

            Spectrum L = emission(light, -dir_light, 0.0, light_sample, scene);
            Spectrum f = eval(*mat, dir_view, dir_light, *pathvertex, scene.texture_pool);
            
            Real pdf_pol = pdf_point_on_light(light, light_sample, p_initial, scene);
            
            Real pdf_nee = light_pmf(scene, sampled_light) * pdf_pol * avg(p_trans_nee);
            Spectrum contrib = T_light * G * f * L / pdf_nee;  

            Real pdf_bsdf = pdf_sample_bsdf(*mat, dir_view, dir_light, *pathvertex, scene.texture_pool) * G * avg(p_trans_dir);
            Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_bsdf * pdf_bsdf);

            return w * contrib;
        } else {


            PhaseFunction phase = get_phase_function(scene.media[current_medium]);
            Real G = max(-dot(dir_light, light_sample.normal), 0.0) / (dist_light * dist_light);

            Spectrum L = emission(light, -dir_light, 0.0, light_sample, scene);
            Spectrum f = eval(phase, dir_view, dir_light);

            Real pdf_pol = pdf_point_on_light(light, light_sample, p_initial, scene);

            Real pdf_nee = light_pmf(scene, sampled_light) * pdf_pol * avg(p_trans_nee);
            Spectrum contrib = T_light * G * f * L / pdf_nee;

            Real pdf_phase = pdf_sample_phase(phase, dir_view, dir_light) * G * avg(p_trans_dir);
            Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
            
            return w * contrib;
        }
    }

    return make_zero_spectrum();    
}


// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {

    // screen width and height
    int w = scene.camera.width;
    int h = scene.camera.height;

    // position in screenspace
    float screenspace_x = (x + next_pcg32_real<Real>(rng)) / w;
    float screenspace_y = (y + next_pcg32_real<Real>(rng)) / h; 
    Vector2 screenspace_pos(screenspace_x, screenspace_y);

    // sample a ray
    Ray ray = sample_primary(scene.camera, screenspace_pos); 

    // "For this homework, we will disable ray diﬀerentials by setting it to 0 and do not change ii."
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)}; 

    int current_medium = scene.camera.medium_id;

    Spectrum radiance = make_zero_spectrum();   

    Spectrum current_path_throughput = make_const_spectrum(1);
    Spectrum multi_trans_pdf = make_const_spectrum(1);
    Spectrum multi_nee_pdf = make_const_spectrum(1);

    Vector3 nee_p_cache{0, 0, 0};

    Real bounces = 0;
    Real dir_pdf = 0; 

    bool never_scatter = true;
    bool never_bsdf = true;


    while (true) {
        bool scatter = false;
        std::optional<PathVertex> intersection = intersect(scene, ray, ray_diff);
        // isect might not intersect a surface, but we might be in a volume
        Real t_hit = INFINITY;

        if (intersection) {
            t_hit = distance(intersection->position, ray.org);
        }
        
        Spectrum transmittance = make_const_spectrum(1);
        Spectrum trans_dir_pdf = make_const_spectrum(1);
        Spectrum trans_nee_pdf = make_const_spectrum(1);
        
        if (current_medium != -1) {
            // sample t s.t. p(t) ~ exp(-sigma_t * t)
            // compute transmittance and trans_pdf
            // if t < t_hit, set scatter = True
            // ...

            const Medium &medium = scene.media[current_medium];
            Spectrum majorant    = get_majorant(medium, ray);

            // sample a channel for sampling
            Real u        = next_pcg32_real<Real>(rng);
            int channel   = std::clamp(int(u * 3), 0, 2);
            Real accum_t  = 0;
            int iteration = 0;

            while (true) {
                if (majorant[channel] <= 0) {
                    break;
                }
                if (iteration >= scene.options.max_null_collisions) {
                    break;
                }

                Real t  = -log(1 - next_pcg32_real<Real>(rng)) / majorant[channel];
                Real dt = t_hit - accum_t;

                // update accumlate distance
                accum_t = min(accum_t + t, t_hit);

                if (t < dt) { // haven't reached the surface
                    // sample from real/fake particle events
                    Vector3 p          = ray.org + accum_t * ray.dir;
                    Spectrum sigma_t   = get_sigma_a(medium, p) + get_sigma_s(medium, p);
                    Spectrum real_prob = sigma_t / majorant;
                    if (next_pcg32_real<Real>(rng) < real_prob[channel]) {
                        // hit a "real" particle
                        scatter       = true;
                        never_scatter = false;

                        transmittance *= exp(-majorant * t) / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * real_prob / max(majorant);
                        
                        ray.org = ray.org + accum_t * ray.dir;
                        // don’t need to account for trans_nee_pdf since we scatter
                        break;
                    } else {
                        // hit a "fake" particle
                        transmittance *= exp(-majorant * t) * (majorant - sigma_t) / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * (1 - real_prob) / max(majorant);
                        trans_nee_pdf *= exp(-majorant * t) * majorant / max(majorant);
                    }
                } else { // reach the surface
                    transmittance *= exp(-majorant * dt);
                    trans_dir_pdf *= exp(-majorant * dt);
                    trans_nee_pdf *= exp(-majorant * dt);
                    ray.org += t_hit * ray.dir;
                    break;
                }
                iteration++;
            }

            multi_trans_pdf *= trans_dir_pdf;
            multi_nee_pdf   *= trans_nee_pdf;

            if (!scatter && !intersection) break;
        } else if (intersection) {
            Real eps = get_intersection_epsilon(scene);
            ray.org  = intersection->position + ray.dir * eps;
        } else {
            break;
        }

        current_path_throughput *= (transmittance / avg(trans_dir_pdf));


        // If we reach a surface and didn’t scatter, include the emission.
        if (!scatter) {
            // reach a surface, include emission 
            if (intersection && is_light(scene.shapes[intersection->shape_id])) {
                if (never_scatter && never_bsdf) {
                    // This is the only way we can see the light source, so
                    // we don’t need multiple importance sampling.
                    radiance += current_path_throughput * emission(*intersection, -ray.dir, scene); 
                } else {
                    // Need to account for next event estimation
                    int light_id = get_area_light_id(scene.shapes[intersection->shape_id]);
                    const Light &light = scene.lights[light_id];
                    Vector3 light_n    = intersection->geometric_normal;
                    Vector3 light_pos  = intersection->position;
                    PointAndNormal light_pn = {light_pos, light_n};

                    // Note that pdf_nee needs to account for the path vertex that issued
                    // next event estimation potentially many bounces ago.
                    // The vertex position is stored in nee_p_cache.
                    // The PDF for sampling the light source using phase function sampling + transmittance sampling
                    // The directional sampling pdf was cached in dir_pdf in solid angle measure.
                    // The transmittance sampling pdf was cached in multi_trans_pdf.

                    Vector3 dir_light = normalize(light_pos - nee_p_cache);
                    Real pdf_nee      = light_pmf(scene, light_id) * pdf_point_on_light(light, light_pn, nee_p_cache, scene) * avg(multi_nee_pdf);
                
                    Real G            = max(-dot(dir_light, intersection->geometric_normal), 0.0) /  distance_squared(nee_p_cache, light_pos);
                    Real real_dir_pdf = dir_pdf * avg(multi_trans_pdf) * G;
                    
                    Real w = (real_dir_pdf * real_dir_pdf) / (real_dir_pdf * real_dir_pdf + pdf_nee * pdf_nee);
                    
                    // current_path_throughput already accounts for transmittance.
                    radiance += current_path_throughput * emission(*intersection, -ray.dir, scene) * w;
                }
                
            }
        } 

        if (bounces == scene.options.max_depth - 1 && scene.options.max_depth != -1) break;

        if (!scatter && intersection) {
            if (intersection->material_id == -1) {
                // index-matching interface, skip through it
                if (intersection->exterior_medium_id != intersection->interior_medium_id) {
                    current_medium = update_medium(*intersection, ray, current_medium);
                }

                Vector3 N = dot(ray.dir, intersection->geometric_normal) > 0 ? intersection->geometric_normal : -intersection->geometric_normal;

                Real eps = get_intersection_epsilon(scene);
                ray.org  = intersection->position + N * eps;

                bounces++;
                continue;
            }
        }

        // sample next direct & update path throughput
        if (scatter){

            never_scatter = false; 

            const Medium &medium = scene.media[current_medium];
            Spectrum sigma_s = get_sigma_s(medium, ray.org);
            Spectrum nee     = next_event_estimation(scene, ray.org, -ray.dir, current_medium, nullptr, nullptr, bounces, rng);
            radiance        += current_path_throughput * sigma_s * nee;

            PhaseFunction phase_func = get_phase_function(medium);
            Vector2 rnd_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir_ = sample_phase_function(phase_func, -ray.dir, rnd_uv);
            current_path_throughput *= eval(phase_func, -ray.dir, *next_dir_) /
                                       pdf_sample_phase(phase_func, -ray.dir, *next_dir_) * sigma_s;

            dir_pdf         = pdf_sample_phase(phase_func, -ray.dir, *next_dir_);
            nee_p_cache     = ray.org;
            multi_trans_pdf = make_const_spectrum(1);
            multi_nee_pdf   = make_const_spectrum(1);

            ray.dir = *next_dir_;
        } else {
            // Hit a surface
            never_bsdf = false;
            if (!intersection.has_value()) break;

            const Material &mat = scene.materials[intersection->material_id];

            Spectrum nee = next_event_estimation(scene, ray.org, -ray.dir, current_medium, &mat, &*intersection, bounces, rng);
            radiance    += current_path_throughput * nee;

            Vector2 bsdf_rnd_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real bsdf_rnd_f = next_pcg32_real<Real>(rng);

            auto bsdf_sample = sample_bsdf(mat, -ray.dir, *intersection, scene.texture_pool, bsdf_rnd_uv, bsdf_rnd_f);
            if (!bsdf_sample) break;

            Vector3 bsdf_dir = bsdf_sample->dir_out;
            Real bsdf_pdf    = pdf_sample_bsdf(mat, -ray.dir, bsdf_dir, *intersection, scene.texture_pool);

            if (bsdf_pdf <= 0) break;

            Spectrum bsdf_f = eval(mat, -ray.dir, bsdf_dir, *intersection, scene.texture_pool);

            current_path_throughput *= bsdf_f / bsdf_pdf;

            if (bsdf_sample->eta != 0) {
                // refraction
                if (intersection->interior_medium_id != intersection->exterior_medium_id) {
                    current_medium = update_medium(*intersection, ray, current_medium);
                }
            }

            nee_p_cache     = ray.org;
            dir_pdf         = bsdf_pdf;
            multi_trans_pdf = make_const_spectrum(1);
            multi_nee_pdf   = make_const_spectrum(1);

            ray.dir = bsdf_dir;

            Vector3 N = dot(ray.dir, intersection->geometric_normal) > 0 ? intersection->geometric_normal : -intersection->geometric_normal;

            Real eps = get_intersection_epsilon(scene);
            ray.org  = intersection->position + N * eps;
        }

        // Russian roulette
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(luminance(current_path_throughput), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;
    }
    return radiance;
}

    