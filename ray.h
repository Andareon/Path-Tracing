#ifndef RAY_TRACING_RAY_H
#define RAY_TRACING_RAY_H

#include "glm/geometric.hpp"
#include "glm/vec2.hpp"
#include "glm/vec3.hpp"
#include "glm/vec4.hpp"

#include "config.h"

class Ray {
private:
    glm::vec4 begin_ = glm::vec4(0);
    glm::vec4 direction_ = glm::vec4(0);
    int depth_ = 0;
    glm::ivec2 coords_ = glm::vec2(0);
    glm::vec3 color_ = glm::vec3(1);

public:
    Ray() = default;
    Ray(glm::vec4 begin, glm::vec4 direction, int depth, glm::ivec2 coords)
            : begin_(begin),
              direction_(normalize(direction)),
              depth_(depth),
              coords_(coords) {}

    glm::vec4 GetBegin() const { return begin_; }

    void SetBegin(glm::vec4 begin) { begin_ = begin; }

    glm::vec4 GetDirection() const { return normalize(direction_); }

    void SetDirection(glm::vec4 direction) { direction_ = direction; }

    int GetDepth() const { return depth_; }

    void SetDepth(int depth) { depth_ = depth; }

    glm::ivec2 GetCoords() const { return coords_; }

    glm::vec3 GetColor() const { return color_; }

    void SetColor(glm::vec3 color) { color_ = color; }

    void Reflect(glm::vec4 begin, glm::vec4 direction, glm::vec3 color) {
        begin_ = begin;
        direction_ = normalize(direction);
        color_ *= color;
        depth_++;
    }

    bool IsValid() const {
        return depth_ < Config::get().max_ray_reflections && color_ != glm::vec3(0);
    }

    void MakeInvalid() { depth_ = Config::get().max_ray_reflections; }
};

#endif