#ifndef RAY_TRACING_CONFIG_H
#define RAY_TRACING_CONFIG_H

#include <string>

class Config {
private:
    Config() = default;

public:
    Config(const Config&) = delete;
    Config& operator=(Config&) = delete;

    int height = 512;
    int width = 512;
    int rays_per_pixel = 20;
    int max_ray_reflections = 8;
    int median = 0;
    int gauss = 0;
    float eps = 1e-6;
    float error = 0.001;
    int update = 32;
    float gamma_correction = 1 / 2.2f;
    std::string model_path = "../models/";
    std::string model_name = "Tor.obj";
    std::string skybox = "";
    static Config& get() {
        static Config instance;
        return instance;
    }

    void set_config(int argc, char* argv[]) {
        for (int i = 1; i < argc - 1; i += 2) {
            if ((std::string)argv[i] == "--H") {
                height = std::atoi(argv[i + 1]);
            }

            if ((std::string)argv[i] == "--W") {
                width = std::atoi(argv[i + 1]);
            }

            if ((std::string)argv[i] == "-RPP") {
                rays_per_pixel = std::atoi(argv[i + 1]);
            }

            if ((std::string)argv[i] == "-MRR") {
                max_ray_reflections = std::atoi(argv[i + 1]);
            }

            if ((std::string)argv[i] == "-EPS") {
                eps = static_cast<float>(std::atof(argv[i + 1]));
            }

            if ((std::string)argv[i] == "-ERR") {
                error = static_cast<float>(std::atof(argv[i + 1]));
            }

            if ((std::string)argv[i] == "-MEDIAN") {
                median = std::atoi(argv[i + 1]);
                gauss = 0;
            }

            if ((std::string)argv[i] == "-UPDATE") {
                update = std::atoi(argv[i + 1]);
            }

            if ((std::string)argv[i] == "-MODEL_PATH") {
                model_path = argv[i + 1];
            }

            if ((std::string)argv[i] == "-MODEL_NAME") {
                model_name = argv[i + 1];
            }

            if ((std::string)argv[i] == "-GAUSS") {
                gauss = std::atoi(argv[i + 1]);
                median = 0;
            }

            if ((std::string)argv[i] == "-GAMMA") {
                gamma_correction = static_cast<float>(std::atof(argv[i + 1]));
            }

            if ((std::string)argv[i] == "-SKYBOX") {
                skybox = argv[i + 1];
            }
        }
    }
};

#endif