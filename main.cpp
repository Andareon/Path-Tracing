#include <deque>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <ctime>
#include <string>
#include <assert.h>
#include <random>

#include "bitmap_image.hpp"
#include "glm/geometric.hpp"
#include "glm/vec3.hpp"

using namespace std;
using namespace glm;

static const vec3 white(1, 1, 1);
static const vec3 red(1, 0, 0);
static const vec3 blue(0, 0, 1);
static const vec3 violet(1, 0, 1);
static const vec3 yellow(1, 1, 0);
static const vec3 black(0, 0, 0);
static const vec3 green(0, 1, 0);

const int W = 500;
const int H = 500;
const int RAYS_PER_PIXEL = 4;
const int MAX_RAY_REFLECTIONS = 4;
const float EPS = 1e-4;

float saturate(float x) {
    return std::max(std::min(1.0f, x), 0.0f);
}

vec3 saturate(vec3 x) {
    return vec3(saturate(x.r), saturate(x.g), saturate(x.b));
}

class Ray {
private:
    vec3 begin;
    vec3 dir;
    int depth = 0;
    ivec2 coords = {0, 0};
    vec3 col = vec3(1);

public:

    vec3 getBegin() const {return begin;}

    void setBegin(vec3 bg) {begin=bg;}

    vec3 getDir() const {return dir;}

    void setDir(vec3 dr) {dir=dr;}

    int getDepth() const {return depth;}

    void setDepth(int dp) {depth=dp;}

    ivec2 getCoords() const {return coords;}
    
    vec3 getCol() const {return col;}

    void setCol(vec3 cl) {col=cl;}

    void reflect(vec3 bg, vec3 dr, vec3 cl) {
        begin = bg;
        dir = dr;
        col.r *= cl.r;
        col.g *= cl.g;
        col.b *= cl.b;
        ++depth;
    }

    bool is_valid() const {
        return depth < MAX_RAY_REFLECTIONS;
    }

    void make_invalid() {
        depth = MAX_RAY_REFLECTIONS + 1;
    }

    Ray(vec3 i, vec3 j, int k, ivec2 l) :begin(i), dir(normalize(j)), depth(k), coords(l) {}
};

class BaseMaterial {
public:
    virtual void process(Ray &ray, vec3 pi, vec3 N, vec3 L, float t, vec3 col, vec3 lc, vector<vector<vec3> > &ColorMap, vector<vector<int> > &samplesCount, deque<Ray> &Rays) = 0;
    virtual ~BaseMaterial() = default;
};

class MirrorMaterial : public BaseMaterial {
public:
    void process(Ray &ray, vec3 pi, vec3 N, vec3 L, float t, vec3 col, vec3 lc, vector<vector<vec3> > &ColorMap, vector<vector<int> > &samplesCount, deque<Ray> &Rays) {
        if (ray.getCol() == black) {
            ray.make_invalid();
            return;
        }

        ray.reflect(pi + N * EPS, reflect(ray.getDir(), N), white);

    }
};

class DiffuseMaterial : public BaseMaterial {
public:
    void process(Ray &ray, vec3 pi, vec3 N, vec3 L, float t, vec3 col, vec3 lc, vector<vector<vec3> > &ColorMap, vector<vector<int> > &samplesCount, deque<Ray> &Rays) {
        if (ray.getCol() == black) {
            ray.make_invalid();
            return;
        }

        std::mt19937 gen(time(0));
        std::uniform_real_distribution<> urd(-100.0f, 100.0f);

        vec3 rnd = {urd(gen), urd(gen), urd(gen)};
        rnd = rnd;
        if (dot(N, rnd) < 0) {
            rnd *= -1;
        }
        ray.reflect(pi + N * EPS, rnd, col);

//        float dt = std::max(dot(normalize(L), normalize(N)), 0.0f);
//        dt += 0.05f;
//        vec3 result = col * dt;
//        ColorMap[ray.getCoords().x][ray.getCoords().y] += saturate(result);
//        samplesCount[ray.getCoords().x][ray.getCoords().y] += 1;


    }
};

class LightMaterial : public BaseMaterial {
public:
    void process(Ray &ray, vec3 pi, vec3 N, vec3 L, float t, vec3 col, vec3 lc, vector<vector<vec3> > &ColorMap, vector<vector<int> > &samplesCount, deque<Ray> &Rays) {
//        vec3 result = multiplyColor(lc, ray.getCol());
//        ColorMap[ray.getCoords().x][ray.getCoords().y] += result;
//        ++samplesCount[ray.getCoords().x][ray.getCoords().y];

        ray.reflect(pi + N * EPS, vec3(0), lc);
        ColorMap[ray.getCoords().x][ray.getCoords().y] += ray.getCol();
        ++samplesCount[ray.getCoords().x][ray.getCoords().y];
        ray.make_invalid();

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

    bool intersect(Ray &ray, float &t) const {
        vec3 o = ray.getBegin();
        vec3 d = ray.getDir();
        vec3 oc = o-c;
        float b = dot(oc, d);
        float c = dot(oc, oc) - r*r;
        float disc = b*b-c;
        if (disc < 0) {
            return false;
        } else {
            disc = sqrt(disc);
            float t0 = -b-disc;
            float t1 = -b+disc;
            assert(t0 * t1 > 0);
            float newT = (t0 > 0) ? t0 :t1;
            if (t > newT && newT > 0) {
                t = newT;
                return true;
            } else {
                return false;
            }
        }
    }
};

struct Plane {
    float A, B, C, D;
    vec3 col;
    vec3 N;
    unique_ptr<BaseMaterial> material;
    Plane(float i, float j, float g, float f, vec3 k, unique_ptr<BaseMaterial> &&m){A=i, B=j, C=g, D=f, col=k,
                                                                                    material=move(m),
                                                                                    N=normalize(vec3(A, B, C));}

    vec3 getNormal(vec3 pi) const {
        return N;
    }

    bool intersect(Ray &ray, float &t) const {
        vec3 o = ray.getBegin();
        vec3 d = ray.getDir();
        if (dot(N, d) == 0) {
            return false;
        } else {
            float newT = -((dot(o, N) + D) / dot(d, N));
            if (t > newT && newT > 0) {
                t = newT;
                return true;
            } else {
                return false;
            }
        }
    }
};

void traceRay(Ray &ray, vector<vector<vec3> > &ColorMap, vector<vector<int> > &samplesCount, deque<Ray> &Rays) {


    const int spheres_count = 4;
    const int planes_count = 1;
    static const Sphere spheres[spheres_count] = {Sphere(vec3(-11,7,20),5, red, make_unique<DiffuseMaterial>()),
                                      Sphere(vec3(0,0,35),5,blue, make_unique<MirrorMaterial>()),
                                      Sphere(vec3(11,0,25),5, violet, make_unique<DiffuseMaterial>()),
                                      Sphere(vec3(0,15,5),2, white, make_unique<LightMaterial>())};

    static const Plane planes[planes_count] = {Plane(0, 0, -1, 40, green, make_unique<DiffuseMaterial>())};

    float t = INFINITY;
    int i = 0, cur=-1, type=0;
    for (i = 0; i < spheres_count; ++i) {
        if (spheres[i].intersect(ray, t)) {
            cur = i;
            type = 1;
        }
    }
    for (i = 0; i < planes_count; ++i) {
        if (planes[i].intersect(ray, t)) {
            cur = i;
            type = 2;
        }
    }

    if (type == 1) {
        vec3 pi = ray.getBegin() + ray.getDir() * t;
        vec3 N = spheres[cur].getNormal(pi);
        vec3 L = spheres[3].c - pi;
        spheres[cur].material->process(ray, pi, N, L, t, spheres[cur].col, spheres[3].col, ColorMap, samplesCount, Rays);
    } else if (type == 2) {
        vec3 pi = ray.getBegin() + ray.getDir() * t;
        vec3 N = planes[cur].getNormal(pi);
        vec3 L = spheres[3].c - pi;
        planes[cur].material->process(ray, pi, N, L, t, planes[cur].col, spheres[3].col, ColorMap, samplesCount, Rays);
    }
}



int main() {
    std::mt19937 gen(time(0));
    std::uniform_real_distribution<> urd(-0.5f, 0.5f);


    unsigned int start_time = clock();


    vector<vector<vec3> > ColorMap(W, vector<vec3>(H, vec3(0)));
    vector<vector<int> > samplesCount(W, vector<int>(H, 0));

    bitmap_image image(H, W);

    deque<Ray> Rays;

    image.clear();


    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            if (x % 100 == 0)
            for (int i = 0; i < RAYS_PER_PIXEL; ++i) {
                vec3 dir = vec3((x - W / 2 + urd(gen)) / W,
                                -(y - H / 2 + urd(gen)) / H,
                                1);
                Ray cur = Ray(vec3(0, 0, -20), dir, 0, ivec2(x, y));
                while (cur.is_valid()) {
                    traceRay(cur, ColorMap, samplesCount, Rays);
                }
            }
        }
    }

    for (int x = 0; x < W; ++x) {
        for (int y = 0; y < H; ++y) {
            vec3 c = ColorMap[x][y] / static_cast<float>(samplesCount[x][y]) * 255.0f;
            image.set_pixel(x, y, c.r, c.g, c.b);
        }
    }

    unsigned int end_time = clock();

    time_t t = time(0);
    struct tm * now = localtime( & t );
    string date = to_string(now->tm_year + 1900) + '-' + to_string(now->tm_mon + 1) + '-' +
                  to_string(now->tm_mday) + '-' + to_string(now->tm_hour) + '-' + to_string(now->tm_min) + '-' +
                  to_string(now->tm_sec) + "  " + to_string(end_time - start_time);
    image.save_image(date + ".bmp");
    image.save_image("last.bmp");
}