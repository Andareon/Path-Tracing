#ifndef RAY_TRACING_MATERIAL_H
#define RAY_TRACING_MATERIAL_H

#include <functional>
#include <random>

#include <ctime>

#include "ray.h"
#include "timer.h"

const float pi = 3.141593;

using MaterialProcessor = std::function<void(Ray &, glm::vec4, glm::vec4)>;

inline float Random() {
    static std::default_random_engine generator(Config::get().getSeed());
    static std::uniform_real_distribution<float> distribution(0.0, 1.0);
    return distribution(generator);
}

struct MaterialCharacteristics {
    glm::vec3 Ke = glm::vec3(0);
    glm::vec3 Kd = glm::vec3(0);
};

class Material {
private:
    std::vector<float> chance_;
    std::vector<MaterialProcessor> functions_;

public:
    void Process(Ray &ray, const glm::vec4& drop_point, const glm::vec4& normal) {
        if (functions_.empty()) {
            ray.MakeInvalid();
        } else if (functions_.size() == 1) {
            functions_[0](ray, drop_point, normal);
        } else {
            float sample = Random();
            int i = 0;
            while (sample > 0) {
                sample -= chance_[i];
                ++i;
            }
            functions_[i](ray, drop_point, normal);
        }
    }

    void AddFunctions(MaterialProcessor function, float chance) {
        functions_.emplace_back(function);
        chance_.emplace_back(chance);
    }
};

inline Material Factory(const MaterialCharacteristics &characteristics,
                        std::vector<std::vector<glm::vec3> > &color_map,
                        std::vector<std::vector<glm::vec3> > &color2_map,
                        std::vector<std::vector<int> > &samples_count) {
    const glm::vec3 Kd = characteristics.Kd;
    const glm::vec3 Ke = characteristics.Ke;
    Material material;
    if (Ke != glm::vec3(0)) {
        auto function = [&color_map, &color2_map, &samples_count, Kd](
                Ray &ray, glm::vec4, glm::vec4 N) {
            if (dot(ray.GetDirection(), N) > 0) {
                ray.MakeInvalid();
                return;
            }
            color_map[ray.GetCoords().x][ray.GetCoords().y] += ray.GetColor() * Kd;
            color2_map[ray.GetCoords().x][ray.GetCoords().y] +=
                    (ray.GetColor() * Kd) * (ray.GetColor() * Kd);
            ++samples_count[ray.GetCoords().x][ray.GetCoords().y];
            ray.MakeInvalid();
        };
        material.AddFunctions(function, 1);
    } else {
        auto function = [Kd](Ray &ray, glm::vec4 drop_point, glm::vec4 N) {
            const float xi1 = Random(), xi2 = Random();
            glm::vec4 rnd =
                    glm::normalize(glm::vec4(std::sqrt(xi1) * std::cos(2 * pi * xi2),
                                             std::sqrt(xi1) * std::sin(2 * pi * xi2), std::sqrt(1 - xi1), 0));
            if (glm::dot(N, rnd) < 0) {
                rnd *= -1;
            }
            const float dt = std::max(0.0f, glm::dot(N, rnd));
            ray.Reflect(drop_point + N * Config::get().eps, rnd, Kd * dt);
        };
        material.AddFunctions(function, 1);
    }
    return material;
}

#endif  // RAY_TRACING_MATERIAL_H