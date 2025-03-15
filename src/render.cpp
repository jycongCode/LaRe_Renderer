#include "render.h"
#include "intersection.h"
#include "material.h"
#include "parallel.h"
#include "path_tracing.h"
#include "vol_path_tracing.h"
#include "restir_path_tracing.h"
#include "pcg.h"
#include "progress_reporter.h"
#include "scene.h"
/// Render auxiliary buffers e.g., depth.
Image3 aux_render(const Scene &scene) {
    int w = scene.camera.width, h = scene.camera.height;
    Image3 img(w, h);

    constexpr int tile_size = 16;
    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;

    parallel_for([&](const Vector2i &tile) {
        int x0 = tile[0] * tile_size;
        int x1 = min(x0 + tile_size, w);
        int y0 = tile[1] * tile_size;
        int y1 = min(y0 + tile_size, h);
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                Ray ray = sample_primary(scene.camera, Vector2((x + Real(0.5)) / w, (y + Real(0.5)) / h));
                RayDifferential ray_diff = init_ray_differential(w, h);
                if (std::optional<PathVertex> vertex = intersect(scene, ray, ray_diff)) {
                    Real dist = distance(vertex->position, ray.org);
                    Vector3 color{0, 0, 0};
                    if (scene.options.integrator == Integrator::Depth) {
                        color = Vector3{dist, dist, dist};
                    } else if (scene.options.integrator == Integrator::ShadingNormal) {
                        // color = (vertex->shading_frame.n + Vector3{1, 1, 1}) / Real(2);
                        color = vertex->shading_frame.n;
                    } else if (scene.options.integrator == Integrator::MeanCurvature) {
                        Real kappa = vertex->mean_curvature;
                        color = Vector3{kappa, kappa, kappa};
                    } else if (scene.options.integrator == Integrator::RayDifferential) {
                        color = Vector3{ray_diff.radius, ray_diff.spread, Real(0)};
                    } else if (scene.options.integrator == Integrator::MipmapLevel) {
                        const Material &mat = scene.materials[vertex->material_id];
                        const TextureSpectrum &texture = get_texture(mat);
                        auto *t = std::get_if<ImageTexture<Spectrum>>(&texture);
                        if (t != nullptr) {
                            const Mipmap3 &mipmap = get_img3(scene.texture_pool, t->texture_id);
                            Vector2 uv{modulo(vertex->uv[0] * t->uscale, Real(1)),
                                       modulo(vertex->uv[1] * t->vscale, Real(1))};
                            // ray_diff.radius stores approximatedly dpdx,
                            // but we want dudx -- we get it through
                            // dpdx / dpdu
                            Real footprint = vertex->uv_screen_size;
                            Real scaled_footprint = max(get_width(mipmap), get_height(mipmap)) *
                                                    max(t->uscale, t->vscale) * footprint;
                            Real level = log2(max(scaled_footprint, Real(1e-8f)));
                            color = Vector3{level, level, level};
                        }
                    }
                    img(x, y) = color;
                } else {
                    img(x, y) = Vector3{0, 0, 0};
                }
            }
        }
    }, Vector2i(num_tiles_x, num_tiles_y));

    return img;
}

Image3 path_render(const Scene &scene) {
    int w = scene.camera.width, h = scene.camera.height;
    Image3 img(w, h);

    constexpr int tile_size = 16;
    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;

    ProgressReporter reporter(num_tiles_x * num_tiles_y);
    parallel_for([&](const Vector2i &tile) {
        // Use a different rng stream for each thread.
        pcg32_state rng = init_pcg32(tile[1] * num_tiles_x + tile[0]);
        int x0 = tile[0] * tile_size;
        int x1 = min(x0 + tile_size, w);
        int y0 = tile[1] * tile_size;
        int y1 = min(y0 + tile_size, h);
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                Spectrum radiance = make_zero_spectrum();
                int spp = scene.options.samples_per_pixel;
                for (int s = 0; s < spp; s++) {
                    radiance += path_tracing(scene, x, y, rng);
                }
                img(x, y) = radiance / Real(spp);
            }
        }
        reporter.update(1);
    }, Vector2i(num_tiles_x, num_tiles_y));
    reporter.done();
    return img;
}

Image3 vol_path_render(const Scene &scene) {
    int w = scene.camera.width, h = scene.camera.height;
    Image3 img(w, h);

    constexpr int tile_size = 16;
    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;

    auto f = vol_path_tracing;
    if (scene.options.vol_path_version == 1) {
        f = vol_path_tracing_1;
    } else if (scene.options.vol_path_version == 2) {
        f = vol_path_tracing_2;
    } else if (scene.options.vol_path_version == 3) {
        f = vol_path_tracing_3;
    } else if (scene.options.vol_path_version == 4) {
        f = vol_path_tracing_4;
    } else if (scene.options.vol_path_version == 5) {
        f = vol_path_tracing_5;
    } else if (scene.options.vol_path_version == 6) {
        f = vol_path_tracing;
    }

    ProgressReporter reporter(num_tiles_x * num_tiles_y);
    parallel_for([&](const Vector2i &tile) {
        // Use a different rng stream for each thread.
        pcg32_state rng = init_pcg32(tile[1] * num_tiles_x + tile[0]);
        int x0 = tile[0] * tile_size;
        int x1 = min(x0 + tile_size, w);
        int y0 = tile[1] * tile_size;
        int y1 = min(y0 + tile_size, h);
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                Spectrum radiance = make_zero_spectrum();
                int spp = scene.options.samples_per_pixel;
                for (int s = 0; s < spp; s++) {
                    Spectrum L = f(scene, x, y, rng);
                    if (isfinite(L)) {
                        // Hacky: exclude NaNs in the rendering.
                        radiance += L;
                    }
                }
                img(x, y) = radiance / Real(spp);
            }
        }
        reporter.update(1);
    }, Vector2i(num_tiles_x, num_tiles_y));
    reporter.done();
    return img;
}

Image3 restir_render(const Scene &scene) {
    int w = scene.camera.width, h = scene.camera.height;
    Image3 img(w, h);
    

    constexpr int tile_size = 16;
    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;

    Reservoir* buffer_1 = new Reservoir[w * h];
    Reservoir* buffer_2 = new Reservoir[w * h];
    Reservoir* shade_buffer = buffer_2;
    for(int i = 0;i<scene.options.N;++i){
        std::cout << "Initial Sampling ..." << std::endl;
        ProgressReporter reporter_1(num_tiles_x * num_tiles_y);
        parallel_for([&](const Vector2i &tile) {
            // Use a different rng stream for each thread.
            pcg32_state rng = init_pcg32(tile[1] * num_tiles_x + tile[0]);
            int x0 = tile[0] * tile_size;
            int x1 = min(x0 + tile_size, w);
            int y0 = tile[1] * tile_size;
            int y1 = min(y0 + tile_size, h);
            for (int y = y0; y < y1; y++) {
                for (int x = x0; x < x1; x++) {
                    initial_sample(scene,buffer_1,x,y,rng);
                }
            }
            reporter_1.update(1);
        }, Vector2i(num_tiles_x, num_tiles_y));
        reporter_1.done();

        std::cout << "Shadow visible ..." << std::endl;
        ProgressReporter reporter_s(num_tiles_x * num_tiles_y);
        parallel_for([&](const Vector2i &tile) {
            // Use a different rng stream for each thread.
            pcg32_state rng = init_pcg32(tile[1] * num_tiles_x + tile[0]);
            int x0 = tile[0] * tile_size;
            int x1 = min(x0 + tile_size, w);
            int y0 = tile[1] * tile_size;
            int y1 = min(y0 + tile_size, h);
            for (int y = y0; y < y1; y++) {
                for (int x = x0; x < x1; x++) {
                    shadow_pixel(scene,buffer_1,x,y,rng);
                }
            }
            reporter_s.update(1);
        }, Vector2i(num_tiles_x, num_tiles_y));
        reporter_s.done();

        for(int i = 0;i<scene.options.spc_iter_num;++i){
            std::cout << "Spatial Reuse " << i << " ..." <<  std::endl;
            ProgressReporter reporter_2(num_tiles_x * num_tiles_y);
            parallel_for([&](const Vector2i &tile) {
                // Use a different rng stream for each thread.
                pcg32_state rng = init_pcg32(tile[1] * num_tiles_x + tile[0]);
                int x0 = tile[0] * tile_size;
                int x1 = min(x0 + tile_size, w);
                int y0 = tile[1] * tile_size;
                int y1 = min(y0 + tile_size, h);
                for (int y = y0; y < y1; y++) {
                    for (int x = x0; x < x1; x++) {
                        if(i % 2){                        
                            spatial_combine(scene,buffer_2,buffer_1,x,y,rng);
                        }else{
                            spatial_combine(scene,buffer_1,buffer_2,x,y,rng);
                        }
                    }
                }
                reporter_2.update(1);
            }, Vector2i(num_tiles_x, num_tiles_y));
            reporter_2.done();
            if(i % 2)shade_buffer = buffer_1;
            else shade_buffer = buffer_2;
        }

        std::cout << "Shade Pixel ..." << std::endl;
        ProgressReporter reporter_3(num_tiles_x * num_tiles_y);
        parallel_for([&](const Vector2i &tile) {
            // Use a different rng stream for each thread.
            pcg32_state rng = init_pcg32(tile[1] * num_tiles_x + tile[0]);
            int x0 = tile[0] * tile_size;
            int x1 = min(x0 + tile_size, w);
            int y0 = tile[1] * tile_size;
            int y1 = min(y0 + tile_size, h);
            for (int y = y0; y < y1; y++) {
                for (int x = x0; x < x1; x++) {
                    img(x,y) += shade_pixel(scene,shade_buffer,x,y,rng);
                    if(i+1 == scene.options.N)
                        img(x,y) /= Real(scene.options.N);
                        // img(x,y)[0] = pow(img(x,y)[0],Real(1)/2.2);
                        img(x,y)[0] = img(x,y)[0] / (img(x,y)[0]+Real(1));
                        // img(x,y)[1] = pow(img(x,y)[1],Real(1)/2.2);
                        img(x,y)[1] = img(x,y)[1] / (img(x,y)[1]+Real(1));
                        // img(x,y)[2] = pow(img(x,y)[2],Real(1)/2.2);
                        img(x,y)[2] = img(x,y)[2] / (img(x,y)[2]+Real(1));
                }
            }
            reporter_3.update(1);
        }, Vector2i(num_tiles_x, num_tiles_y));
        reporter_3.done();
    }
    return img;
}

Image3 render(const Scene &scene) {
    if (scene.options.integrator == Integrator::Depth ||
            scene.options.integrator == Integrator::ShadingNormal ||
            scene.options.integrator == Integrator::MeanCurvature ||
            scene.options.integrator == Integrator::RayDifferential ||
            scene.options.integrator == Integrator::MipmapLevel) {
        return aux_render(scene);
    } else if (scene.options.integrator == Integrator::Path) {
        return path_render(scene);
    } else if (scene.options.integrator == Integrator::VolPath) {
        return vol_path_render(scene);
    } else if(scene.options.integrator == Integrator::ReSTIR){
        return restir_render(scene);
    } else {
        assert(false);
        return Image3();
    }
}
