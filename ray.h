#ifndef RAY_TRACING_RAY_H
#define RAY_TRACING_RAY_H

#include "glm/geometric.hpp"

using namespace glm;

class Ray {
private:
    vec4 begin;
    vec4 dir;
    int depth = 0;
    ivec2 coords = {0, 0};
    vec3 col = vec3(1);

public:
    vec4 getBegin() const {return begin;}

    void setBegin(vec4 bg) {begin=bg;}

    vec4 getDir() const {return normalize(dir);}

    void setDir(vec4 dr) {dir=dr;}

    int getDepth() const {return depth;}

    void setDepth(int dp) {depth=dp;}

    ivec2 getCoords() const {return coords;}

    vec3 getCol() const {return col;}

    void setCol(vec3 cl) {col=cl;}

    void reflect(vec4 bg, vec4 dr, vec3 cl) {
        begin = bg;
        dir = normalize(dr);
        col *= cl;
        depth++;
    }

    bool is_valid() const {return depth < Config::get().MAX_RAY_REFLECTIONS;}

    void make_invalid() {depth = Config::get().MAX_RAY_REFLECTIONS;}

    Ray(vec4 i, vec4 j, int k, ivec2 l) :begin(i), dir(normalize(j)), depth(k), coords(l) {}
};

#endif