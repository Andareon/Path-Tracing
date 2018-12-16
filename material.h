#ifndef RAY_TRACING_MATERIAL_H
#define RAY_TRACING_MATERIAL_H

#include <functional>
#include <random>

#include "ray.h"

const float PI = 3.141593;

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

Material Factory(std::vector<float> args, std::vector<std::vector<glm::vec3> > &ColorMap, std::vector<std::vector<glm::vec3> > &Color2Map, std::vector<std::vector<int> > &SamplesCount) {
    glm::vec3 Kd = glm::vec3(args[0], args[1], args[2]);
    glm::vec3 Ke = glm::vec3(args[3], args[4], args[5]);
    Material mat;
    if (Ke != glm::vec3(0)) {
        auto func = [&ColorMap, &Color2Map, &SamplesCount, Kd] (Ray &ray, glm::vec4 drop_point, glm::vec4 N) {
            if (dot(ray.getDir(), N) < 0) {
                ray.make_invalid();
                return;
            }
            ColorMap[ray.getCoords().x][ray.getCoords().y] += ray.getCol() * Kd;
            Color2Map[ray.getCoords().x][ray.getCoords().y] += (ray.getCol() * Kd) * (ray.getCol() * Kd);
            ++SamplesCount[ray.getCoords().x][ray.getCoords().y];
            ray.make_invalid();};
        mat.add_functions(func, 1);
    } else {
        auto func = [Kd] (Ray &ray, glm::vec4 drop_point, glm::vec4 N) {
            float xi1 = random0_1(),
                    xi2 = random0_1();
            glm::vec4 rnd = normalize(glm::vec4(sqrt(xi1) * cos(2 * PI * xi2), sqrt(xi1) * sin(2 * PI * xi2), sqrt(1 - xi1), 0));
            if (dot(N, rnd) < 0) {
                rnd *= -1;
            }
            float dt = std::max(0.0f, dot(N, rnd));
            ray.reflect(drop_point + N * Config::get().EPS, rnd, Kd * dt);
        };
        mat.add_functions(func, 1);
    }
    return mat;
}

#endif //RAY_TRACING_MATERIAL_H
