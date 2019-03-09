#include <utility>

#ifndef RAY_TRACING_TRIANGLES_H
#define RAY_TRACING_TRIANGLES_H

#include "glm/geometric.hpp"

#include "material.h"

inline float PlaneIntersect(const Ray &ray, const glm::vec4& plane) {
    const float signedDistToPlane = ray.GetDirection().x * plane.x + ray.GetDirection().y * plane.y + ray.GetDirection().z * plane.z;
    return -(ray.GetBegin().x * plane.x + ray.GetBegin().y * plane.y + ray.GetBegin().z * plane.z + plane.w) / signedDistToPlane;
}

inline float ParallelogramSquare(const glm::vec3 &dir1, const glm::vec3 &dir2) {
    return glm::length(glm::cross(dir1, dir2));
}

class Triangle {
private:
    glm::vec4 plane_;
    Material material_;
    std::array<glm::vec3, 3> vertices_;
    float square;

public:
    Triangle(const std::array<glm::vec4, 3> &vertices, Material material)
            : material_(std::move(material)) {
        for (int i = 0; i < vertices.size(); ++i) {
            vertices_[i] = vertices[i];
        }
        glm::vec3 AB = vertices_[1] - vertices_[0];
        glm::vec3 AC = vertices_[2] - vertices_[0];
        SetNormal(glm::cross(AB, AC));
        square = ParallelogramSquare(AB, AC);
    }

    glm::vec4 GetNormal() const { return glm::vec4(plane_.x, plane_.y, plane_.z, 0); }

    void SetNormal(glm::vec3 normal) {
        normal = normalize(normal);
        plane_ = glm::vec4(normal, 0);
        plane_.w = -dot(normal, vertices_[0]);
    }

    Material GetMaterial() const { return material_; }

    bool Intersect(const Ray &ray, float &distance) const {
        const float new_distance = PlaneIntersect(ray, plane_);

        if (new_distance >= distance || new_distance < Config::get().eps) {
            return false;
        }

        const glm::vec3 drop_point = ray.GetBegin() + ray.GetDirection() * new_distance;
        const glm::vec3 fromVert0 = drop_point - vertices_[0];
        const glm::vec3 fromVert1 = drop_point - vertices_[1];
        const glm::vec3 fromVert2 = drop_point - vertices_[2];
        const float small_square_1 = ParallelogramSquare(fromVert0, fromVert1);
        if (small_square_1 > square + Config::get().eps) {
            return false;
        }
        const float small_square_2 = ParallelogramSquare(fromVert0, fromVert2);
        if (small_square_1 + small_square_2 > square + Config::get().eps) {
            return false;
        }
        const float small_square_3 = ParallelogramSquare(fromVert2, fromVert1);
        if (std::abs(square - small_square_1 - small_square_2 - small_square_3) > Config::get().eps) {
            return false;
        }
        distance = new_distance;
        return true;
    }
};


#endif //RAY_TRACING_TRIANGLES_H