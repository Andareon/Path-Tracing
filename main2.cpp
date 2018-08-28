#include <deque>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "bitmap_image.hpp"
#include "glm/geometric.hpp"
#include "glm/vec3.hpp"

using namespace std;
using namespace glm;

struct Ray;
struct BaseMaterial;
struct MirrorMaterial;
struct Sphere;
float randFloat(const float min, const float max);
float saturate(float x);
vec3 traceRay(const Ray& ray);



static const vec3 white(255, 255, 255);
static const vec3 red(255, 0, 0);
static const vec3 blue(0, 0, 255);
static const vec3 violet(255, 0, 255);
static const vec3 yellow(255, 255, 0);
static const vec3 black(0, 0, 0);

const int W = 500;
const int H = 500;
const int RAYS_PER_PIXEL = 2;
const int MAX_RAY_REFLECTIONS = 4;

deque<Ray> Deq;

struct Ray {
    vec3 o;
    vec3 d;
    int depth = 0;
    int x=0;
    int y=0;
    Ray(vec3 i, vec3 j, int k, int l, int m) {o=i, d=j, depth=k, x=l, y=m;}
};

class BaseMaterial {
public:
    virtual void process(Ray ray, vec3 pi, vec3 N, vec3 L, vec3 &result, float t, vec3 col, vec3 lc) = 0;
    virtual ~BaseMaterial() = default;
};

class MirrorMaterial : public BaseMaterial {
public:
    void process(Ray ray, vec3 pi, vec3 N, vec3 L, vec3 &result, float t, vec3 col, vec3 lc) {
        vec3 ans = reflect(-ray.d, N);
        Deq.push_back(Ray(pi + N * 1e-7f, normalize(ans), ray.depth + 1, ray.x, ray.y));
    }
};

class DiffuseMaterial : public BaseMaterial {
public:
    void process(Ray ray, vec3 pi, vec3 N, vec3 L, vec3 &result, float t, vec3 col, vec3 lc) {
        float dt = std::max(dot(normalize(L), normalize(N)), 0.0f);
        dt += 0.05f;
        result = col * dt;
    }
};

class LightMaterial : public BaseMaterial {
public:
    void process(Ray ray, vec3 pi, vec3 N, vec3 L, vec3 &result, float t, vec3 col, vec3 lc) {
        result = lc;
    }
};

struct Sphere {
    vec3 c;
    float r;
    vec3 col;
    unique_ptr<BaseMaterial> material;
    Sphere(vec3 i, float j, vec3 k, unique_ptr<BaseMaterial> &&m){c=i, r=j, col=k, material=move(m);}

    vec3 getNormal(vec3 pi) const {
        return (pi - c) / r;
    }

    bool intersect(const Ray &ray, float &t) const {
        vec3 o = ray.o;
        vec3 d = normalize(ray.d);
        vec3 oc = o-c;
        float b = 2 * dot(oc, d);
        float c = dot(oc, oc) - r*r;
        float disc = b*b-4*c;
        if (disc < 0) {
            return false;
        } else {
            disc = sqrt(disc);
            float t0 = (-b-disc)/2;
            float t1 = (-b+disc)/2;
            float newT = (t0 < t1) ? t0 :t1;
            if (t > newT) {
                t = newT;
                return true;
            } else {
                return false;
            }
        }
    }
};

float randFloat(const float min, const float max) {
    return (static_cast<float>(rand()) / RAND_MAX) * (max - min) + min;
}

float saturate(float x) {
    return std::max(std::min(255.0f, x), 0.0f);
}

vec3 saturate(vec3 x) {
    return vec3(saturate(x.r), saturate(x.g), saturate(x.b));
}

vec3 traceRay(const Ray &ray) {


    static const Sphere spheres[4] = {Sphere(vec3(-11,7,20),5, red, make_unique<DiffuseMaterial>()),
                                      Sphere(vec3(0,0,35),5,blue, make_unique<MirrorMaterial>()),
                                      Sphere(vec3(11,0,25),5, violet, make_unique<DiffuseMaterial>()),
                                      Sphere(vec3(0,15,15),2, white, make_unique<LightMaterial>())};

    vec3 result = black;

    float t = 20000000;
    int i = 0, cur=-1;
    for (const auto &sphere: spheres) {
        if (sphere.intersect(ray, t)) {
            cur = i;
        }
        ++i;
    }
    if (cur > -1) {
        vec3 pi = ray.o + ray.d * t;
        vec3 N = spheres[cur].getNormal(pi);
        vec3 L = spheres[3].c - pi;
        spheres[cur].material->process(ray, pi, N, L, result, t, spheres[cur].col, spheres[3].col);
    }
    return saturate(result);
}



int main() {
    vector<vector<vec3> > ColorMap(H, vector<vec3>(W, vec3(0)));
    vector<vector<int> > samplesCount(H, vector<int>(W, 0));

    deque<Ray> Deq;

    bitmap_image image(H, W);
    image.clear();
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            float r = 0, g = 0, b = 0;
            for (int i = 0; i < RAYS_PER_PIXEL; ++i) {
                vec3 dir = vec3((x - W / 2 + randFloat(-0.5f, 0.5f)) / W,
                                -(y - H / 2 + randFloat(-0.5f, 0.5f)) / H,
                                1);
                Deq.push_back(Ray(vec3(0, 0, -20), normalize(dir), 0, x, y));
            }
        }
    }

    while (!Deq.empty()) {
        Ray ray = Deq[0];
        Deq.pop_front();
        if (ray.depth < MAX_RAY_REFLECTIONS) {
            vec3 c = traceRay(ray);
            ColorMap[ray.y][ray.x] += c;
            samplesCount[ray.y][ray.x] += 1;
        }
    }

    for (int x = 0; x < W; ++x) {
        for (int y = 0; y < H; ++y) {
            vec3 c = ColorMap[y][x];
            image.set_pixel(x, y, c.r / samplesCount[y][x], c.g / samplesCount[y][x], c.b / samplesCount[y][x]);
        }
    }

    image.save_image("2.bmp");
}