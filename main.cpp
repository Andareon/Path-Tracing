#include <chrono>

#include "bitmap_image.hpp"

#include "material.h"
#include "scene.h"

using namespace std;
using namespace glm;

vector<vector<vec3> > GaussBlur(vector<vector<vec3> > color_map, float r) {
    const int rs = ceil(r * 2.57);
    vector<vector<vec3> > ans(Config::get().width,
                              vector<vec3>(Config::get().height, vec3(0)));
    for (int i = 0; i < Config::get().height; ++i) {
        for (int j = 0; j < Config::get().width; ++j) {
            vec3 val = vec3(0);
            float wsum = 0;
            for (int iy = i - rs; iy <= i + rs; ++iy) {
                for (int ix = j - rs; ix <= j + rs; ++ix) {
                    int x = std::min(Config::get().width - 1, std::max(0, ix));
                    int y = std::min(Config::get().height - 1, std::max(0, iy));
                    int dsq = (ix - j) * (ix - j) + (iy - i) * (iy - i);
                    float wght = exp(-dsq / (2 * r * r)) / (pi * 2 * r * r);
                    val += color_map[x][y] * wght;
                    wsum += wght;
                }
            }
            ans[j][i] = round(val / wsum);
        }
    }
    return ans;
}

vec4 Refract(vec4 incidence_direction, vec4 normal, float refractive) {
    float cos_incidence = -dot(normal, incidence_direction);
    float sin_refraction =
            refractive * sqrt(1.0f - cos_incidence * cos_incidence);
    if (sin_refraction > 1) {
        return reflect(incidence_direction, normal);
    }
    float cos_refraction = sqrt(1.0f - sin_refraction * sin_refraction);
    vec4 refraction_direction =
            incidence_direction * refractive +
            normal * (refractive * cos_incidence - cos_refraction);
    return refraction_direction;
}

vector<vector<vec3> > MedianFilter(vector<vector<vec3> > color_map,
                                   int window_size) {
    vector<vector<vec3> > ans(Config::get().width,
                              vector<vec3>(Config::get().height, vec3(0)));
    for (int y = 0; y < Config::get().height; ++y) {
        for (int x = 0; x < Config::get().width; ++x) {
            vector<float> window_r;
            vector<float> window_g;
            vector<float> window_b;
            for (int window_x = -window_size; window_x < window_size + 1;
                 ++window_x) {
                for (int window_y = -window_size; window_y < window_size + 1;
                     ++window_y) {
                    int i = std::max(std::min(window_x + x, Config::get().width - 1), 0);
                    int j = std::max(std::min(window_y + y, Config::get().height - 1), 0);
                    window_r.push_back(color_map[i][j].r);
                    window_g.push_back(color_map[i][j].g);
                    window_b.push_back(color_map[i][j].b);
                }
            }
            sort(window_r.begin(), window_r.end());
            ans[x][y].r = window_r[window_size * window_size / 2];

            sort(window_g.begin(), window_g.end());
            ans[x][y].g = window_g[window_size * window_size / 2];

            sort(window_b.begin(), window_b.end());
            ans[x][y].b = window_b[window_size * window_size / 2];
        }
    }
    return ans;
}

int main(int argc, char *argv[]) {
    chrono::milliseconds start_time = chrono::duration_cast<chrono::milliseconds>(
            chrono::system_clock::now().time_since_epoch());

    Config::get().set_config(argc, argv);
    static default_random_engine generator(Config::get().getSeed());
    static uniform_real_distribution<> distribution(-0.5f, 0.5f);

    vector<vector<vec3> > color_map(Config::get().width,
                                    vector<vec3>(Config::get().height, vec3(0)));
    vector<vector<vec3> > color2_map(Config::get().width,
                                     vector<vec3>(Config::get().height, vec3(0)));
    vector<vector<int> > samples_count(Config::get().width,
                                       vector<int>(Config::get().height, 0));

    Scene scene = Scene(color_map, color2_map, samples_count);
    scene.LoadModel(Config::get().model_path + Config::get().model_name);

    bitmap_image image(Config::get().width, Config::get().height);

    image.clear();
    std::vector<Ray> rays;
    rays.resize(Config::get().height * Config::get().width);
    for (int i = 0; i < Config::get().rays_per_pixel; ++i) {
#pragma omp parallel for num_threads(4)
        for (int y = 0; y < Config::get().height; ++y) {
            for (int x = 0; x < Config::get().width; ++x) {
                const auto sample_count = static_cast<float>(samples_count[x][y]);
                if (i > 10 && sample_count > 0) {
                    vec3 D = (color2_map[x][y] / sample_count -
                              (color_map[x][y] / sample_count) *
                              (color_map[x][y] / sample_count));
                    if ((i % 4) && D.r < Config::get().error &&
                        D.g < Config::get().error && D.b < Config::get().error) {
                        continue;
                    }
                }
                vec4 direction = vec4(
                        (x + distribution(generator)) / Config::get().width - 0.5f,
                        -(y + distribution(generator)) / Config::get().height + 0.5f, 1, 0);
                rays[y * Config::get().width + x] = Ray(vec4(0, 0, -20, 1), direction, 0, ivec2(x, y));
            }
        }
#pragma omp parallel for num_threads(4)
        for (int y = 0; y < Config::get().height; ++y) {
            for (int x = 0; x < Config::get().width; ++x) {
                Ray &ray = rays[y * Config::get().width + x];
                while (ray.IsValid()) {
                    scene.TraceRay(ray);
                }
            }
        }
#pragma omp parallel for num_threads(4)
        for (int y = 0; y < Config::get().height; ++y) {
            for (int x = 0; x < Config::get().width; ++x) {
                if (i % Config::get().update == 0) {
                    if (samples_count[x][y]) {
                        vec3 c =
                                pow(color_map[x][y] / static_cast<float>(samples_count[x][y]),
                                    vec3(Config::get().gamma_correction)) *
                                255.0f;
                        image.set_pixel(x, y, c.r, c.g, c.b);
                    }
                }
            }
        }
        if (i % Config::get().update == 0) {
            image.save_image("result.bmp");
            cerr << "Image update" << endl;
        }
        cerr << i + 1 << " rays per pixel were sent" << endl;
    }

    for (int y = 0; y < Config::get().height; ++y) {
        for (int x = 0; x < Config::get().width; ++x) {
            if (!samples_count[x][y]) {
                continue;
            }
            color_map[x][y] =
                    pow(color_map[x][y] / static_cast<float>(samples_count[x][y]),
                        vec3(1 / 2.2f)) *
                    255.0f;
        }
    }

    if (Config::get().gauss) {
        color_map = GaussBlur(color_map, Config::get().gauss);
    }
    if (Config::get().median) {
        color_map = MedianFilter(color_map, Config::get().median);
    }
    for (int y = 0; y < Config::get().height; ++y) {
        for (int x = 0; x < Config::get().width; ++x) {
            if (!samples_count[x][y]) {
                continue;
            }
            vec3 color = color_map[x][y];
            image.set_pixel(x, y, color.r, color.g, color.b);
        }
    }

    chrono::milliseconds end_time = chrono::duration_cast<chrono::milliseconds>(
            chrono::system_clock::now().time_since_epoch());
    time_t t = time(nullptr);
    struct tm *now = localtime(&t);
    string date = to_string(now->tm_year + 1900) + '-' +
                  to_string(now->tm_mon + 1) + '-' + to_string(now->tm_mday) +
                  '-' + to_string(now->tm_hour) + '-' + to_string(now->tm_min) +
                  '-' + to_string(now->tm_sec) + "  " +
                  to_string((end_time - start_time).count()) + "   " +
                  to_string(Config::get().rays_per_pixel);
    image.save_image(date + ".bmp");
    image.save_image("result.bmp");
}