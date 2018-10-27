#ifndef RAY_TRACING_CONFIG_H
#define RAY_TRACING_CONFIG_H

using namespace std;

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
    float EPS = 1e-6;
    float error = 0.001;
    static Config& get() {

        static Config instance;
        return instance;
    }


    void set_config(int argc, char* argv[]) {
        for (int i = 1; i < argc - 1; i += 2) {
            if ((string)argv[i] == "-H") {
                height = atoi(argv[i + 1]);
            }

            if ((string)argv[i] == "-W") {
                width = atoi(argv[i + 1]);
            }

            if ((string)argv[i] == "-RPP") {
                RAYS_PER_PIXEL = atoi(argv[i + 1]);
            }

            if ((string)argv[i] == "-MRR") {
                MAX_RAY_REFLECTIONS = atoi(argv[i + 1]);
            }

            if ((string)argv[i] == "-EPS") {
                EPS = static_cast<float>(atof(argv[i + 1]));
            }

            if ((string)argv[i] == "-ERR") {
                error = static_cast<float>(atof(argv[i + 1]));
            }
        }
    }
};

#endif