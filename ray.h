#ifndef RAY_TRACING_RAY_H
#define RAY_TRACING_RAY_H

#include "glm/geometric.hpp"
#include "glm/vec2.hpp"
#include "glm/vec3.hpp"
#include "glm/vec4.hpp"

#include "config.h"

class Ray {
private:
    glm::vec4 begin;
    glm::vec4 dir;
    int depth = 0;
    glm::ivec2 coords = {0, 0};
    glm::vec3 col = glm::vec3(1);

public:
    glm::vec4 getBegin() const {return begin;}

    void setBegin(glm::vec4 bg) {begin=bg;}

    glm::vec4 getDir() const {return normalize(dir);}

    void setDir(glm::vec4 dr) {dir=dr;}

    int getDepth() const {return depth;}

    void setDepth(int dp) {depth=dp;}

    glm::ivec2 getCoords() const {return coords;}

    glm::vec3 getCol() const {return col;}

    void setCol(glm::vec3 cl) {col=cl;}

    void reflect(glm::vec4 bg, glm::vec4 dr, glm::vec3 cl) {
        begin = bg;
        dir = normalize(dr);
        col *= cl;
        depth++;
    }

    bool is_valid() const {return depth < Config::get().MAX_RAY_REFLECTIONS && col != glm::vec3(0);}

    void make_invalid() {depth = Config::get().MAX_RAY_REFLECTIONS;}

    Ray(glm::vec4 i, glm::vec4 j, int k, glm::ivec2 l) :begin(i), dir(normalize(j)), depth(k), coords(l) {}
};

#endif