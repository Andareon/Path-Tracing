//
// Created by alex on 05.06.2019.
//

#ifndef RAY_TRACING_TRACER_H
#define RAY_TRACING_TRACER_H

#include <vector>

#include "triangles.h"

struct IntersectionOptions {
    glm::vec4 intersectPoint_;
    glm::vec4 normal_;
    int materialIndex_;
};

class BaseTracer {
public:
    virtual bool Trace(Ray &ray, IntersectionOptions &options) = 0;
    virtual ~BaseTracer() = default;
};

class SimpleTracer : public BaseTracer {
    const std::vector<Triangle> &triangles;
public:
    SimpleTracer(const std::vector<Triangle> &tri) : triangles(tri) {}

    bool Trace(Ray &ray, IntersectionOptions &options) final {
        float distance = INFINITY;
        int current_triangle = -1;
        for (int i = 0; i < triangles.size(); ++i) {
            if (triangles[i].Intersect(ray, distance)) {
                current_triangle = i;
            }
        }
        if (current_triangle < 0) {
            return false;
        }
        options.intersectPoint_= ray.GetBegin() + ray.GetDirection() * distance;
        options.normal_ = triangles[current_triangle].GetNormal();
        options.materialIndex_ = triangles[current_triangle].GetMaterial();
        return true;
    }
};

#endif //RAY_TRACING_TRACER_H
