#ifndef RAY_TRACING_SCENE_H
#define RAY_TRACING_SCENE_H

#include <memory>

#include "bitmap_image.hpp"

#include "triangles.h"

class BaseTracer;

class Scene {
private:
    std::vector<Triangle> triangles_;
    std::vector<std::vector<glm::vec3> > &color_map_;
    std::vector<std::vector<glm::vec3> > &color2_map_;
    std::vector<std::vector<int> > &samples_count_;
    std::vector<Material> materials_;
    bitmap_image skybox_;
    std::unique_ptr<BaseTracer> tracer;

public:
    std::vector<int> trianglesIndex_;
    Scene(std::vector<std::vector<glm::vec3> > &color_map, std::vector<std::vector<glm::vec3> > &color2_map,
          std::vector<std::vector<int> > &samples_count);
    void LoadModel(const std::string &path);

    void TraceRay(Ray &ray);
};

#endif //RAY_TRACING_SCENE_H