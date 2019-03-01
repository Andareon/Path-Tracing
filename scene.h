#ifndef RAY_TRACING_SCENE_H
#define RAY_TRACING_SCENE_H

#include "triangles.h"
#include "bitmap_image.hpp"

class Scene {
private:
    std::vector<Triangle> triangles_;
    std::vector<std::vector<glm::vec3> > &color_map_;
    std::vector<std::vector<glm::vec3> > &color2_map_;
    std::vector<std::vector<int> > &samples_count_;
    std::vector<Material> materials_;
    bitmap_image skybox_;

public:
    Scene(std::vector<std::vector<glm::vec3> > &color_map, std::vector<std::vector<glm::vec3> > &color2_map,
          std::vector<std::vector<int> > &samples_count);
    void LoadModel(std::string path);

    void AddTriangle(Triangle triangle);

    void TraceRay(Ray &ray);
};

#endif //RAY_TRACING_SCENE_H