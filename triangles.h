#include <utility>

#ifndef RAY_TRACING_TRIANGLES_H
#define RAY_TRACING_TRIANGLES_H

#include "glm/geometric.hpp"

#include "material.h"

inline bool PlaneIntersect(Ray &ray, float &distance, glm::vec4 plane) {
    glm::vec4 begin = ray.GetBegin();
    glm::vec4 direction = ray.GetDirection();
    if (std::abs(glm::dot(plane, direction)) < Config::get().eps) {
        return false;
    } else {
        float new_distance = -glm::dot(begin, plane) / glm::dot(direction, plane);
        if (distance > new_distance && new_distance > 0) {
            distance = new_distance;
            return true;
        } else {
            return false;
        }
    }
}

inline glm::vec4 Cross(glm::vec4 a, glm::vec4 b) {
    return glm::vec4(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                     a.x * b.y - a.y * b.x, 0);
}

inline float Square(glm::vec4 A, glm::vec4 B, glm::vec4 C) {
    glm::vec4 AB = B - A;
    glm::vec4 AC = C - A;
    glm::vec4 cross = Cross(AB, AC);
    return length(cross) / 2;
}

class Triangle {
private:
    glm::vec4 plane_;
    Material material_;
    std::array<glm::vec4, 3> vertices_;

public:
    Triangle(const std::array<glm::vec4, 3> &vertices, Material material)
            : vertices_(vertices), material_(std::move(material)) {
        glm::vec4 AB = vertices_[1] - vertices_[0];
        glm::vec4 AC = vertices_[2] - vertices_[0];
        plane_ = normalize(Cross(AB, AC));
        plane_.w = -dot(plane_, vertices_[0]);
    }

    glm::vec4 GetNormal() const { return glm::vec4(plane_.x, plane_.y, plane_.z, 0); }

    void SetNormal(glm::vec3 normal) {
        normal = normalize(normal);
        plane_.x = normal.x;
        plane_.y = normal.y;
        plane_.z = normal.z;
        plane_.w = 0;
        plane_.w = -dot(plane_, vertices_[0]);
    }

    Material GetMaterial() const { return material_; }

    bool Intersect(Ray &ray, float &distance) const {
        float new_distance = INFINITY;
        if (!PlaneIntersect(ray, new_distance, plane_)) {
            return false;
        }

        if (new_distance >= distance) {
            return false;
        }

        glm::vec4 drop_point = ray.GetBegin() + ray.GetDirection() * new_distance;
        float full_square = Square(vertices_[0], vertices_[1], vertices_[2]);
        float small_square_1 = Square(drop_point, vertices_[1], vertices_[2]);
        float small_square_2 = Square(vertices_[0], drop_point, vertices_[2]);
        float small_square_3 = Square(vertices_[0], vertices_[1], drop_point);
        if (std::abs(full_square - small_square_1 - small_square_2 - small_square_3) >
            Config::get().eps) {
            return false;
        }
        distance = new_distance;
        return true;
    }
};


#endif //RAY_TRACING_TRIANGLES_H