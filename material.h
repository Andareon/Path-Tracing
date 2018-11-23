#ifndef RAY_TRACING_MATERIAL_H
#define RAY_TRACING_MATERIAL_H

#include "glm/geometric.hpp"
#include "ray.h"
#include <vector>
#include <random>
#include <ctime>
#include <functional>

using MaterialProcessor = std::function<void(Ray &, glm::vec4 , glm::vec4)>;

float random0_1() {
    static std::default_random_engine generator(time(0));
    static std::uniform_real_distribution<float> distribution(0.0, 1.0);
    return distribution(generator);
}

class Material {
private:
    std::vector<float> chance;
    std::vector<MaterialProcessor> func;
public:
    void process(Ray &ray, glm::vec4 drop_point, glm::vec4 N) {
        if (!func.size()) {
            ray.make_invalid();
        } else if (func.size() == 1) {
            func[0](ray, drop_point, N);
        } else {
            float sample = random0_1();
            int i = 0;
            while (sample > 0) {
                sample -= chance[i];
                ++i;
            }
            func[i](ray, drop_point, N);
        }
    }

    void add_functions(MaterialProcessor fnc, float ch) {
        func.emplace_back(fnc);
        chance.emplace_back(ch);
    }
};

#endif //RAY_TRACING_MATERIAL_H
