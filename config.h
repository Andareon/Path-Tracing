#ifndef RAY_TRACING_CONFIG_H
#define RAY_TRACING_CONFIG_H

#include <string>

class Config {
private:
    Config(){}
    Config(const Config&) = delete;
    Config& operator=(Config&) = delete;

public:
    int height = 512;
    int width = 512;
    int RAYS_PER_PIXEL = 20;
    int MAX_RAY_REFLECTIONS = 8;
    int median = 0;
    int gauss = 0;
    float EPS = 1e-6;
    float error = 0.001;
    int update = 32;
    float gamma_kor = 1/2.2;
    std::string model_path = "../";
    std::string model_name = "3.obj";
    static Config& get() {

        static Config instance;
        return instance;
    }


    void set_config(int argc, char* argv[]) {
        for (int i = 1; i < argc - 1; i += 2) {
            if ((std::string)argv[i] == "-H") {
                height = std::atoi(argv[i + 1]);
            }

            if ((std::string)argv[i] == "-W") {
                width = std::atoi(argv[i + 1]);
            }

            if ((std::string)argv[i] == "-RPP") {
                RAYS_PER_PIXEL = std::atoi(argv[i + 1]);
            }

            if ((std::string)argv[i] == "-MRR") {
                MAX_RAY_REFLECTIONS = std::atoi(argv[i + 1]);
            }

            if ((std::string)argv[i] == "-EPS") {
                EPS = static_cast<float>(std::atof(argv[i + 1]));
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
                gamma_kor = static_cast<float>(std::atof(argv[i + 1]));
            }
        }
    }
};

#endif