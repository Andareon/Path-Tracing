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
static const vec3 green(0.3, 1, 0.3);

const int W = 500;
const int H = 500;
const int RAYS_PER_PIXEL = 20;
const int MAX_RAY_REFLECTIONS = 4;
const float EPS = 1e-6;
const int update_img = 10;

std::default_random_engine generator(time(0));
std::uniform_int_distribution<int> distribution1(-100, 100);
std::uniform_real_distribution<> distribution2(-0.5f, 0.5f);

class Ray;
void traceRay(Ray &ray, vector<vector<vec3> > &ColorMap, vector<vector<int> > &samplesCount);

float saturate(float x) {
    return std::max(std::min(1.0f, x), 0.0f);
}

vec3 saturate(vec3 x) {
    return vec3(saturate(x.r), saturate(x.g), saturate(x.b));
}

float square(vec3 A, vec3 B, vec3 C) {
    vec3 a = B - A;
    vec3 b = C - A;
    vec3 c = cross(a, b);
    return sqrt(c.x * c.x + c.y * c.y + c.z * c.z) / 2;
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

    vec3 getDir() const {return normalize(dir);}

    void setDir(vec3 dr) {dir=dr;}

    int getDepth() const {return depth;}

    void setDepth(int dp) {depth=dp;}

    ivec2 getCoords() const {return coords;}
    
    vec3 getCol() const {return col;}

    void setCol(vec3 cl) {col=cl;}

    void reflect(vec3 bg, vec3 dr, vec3 cl) {
        begin = bg;
        dir = normalize(dr);
        col *= cl;
        depth++;
    }

    bool is_valid() const {return depth < MAX_RAY_REFLECTIONS;}

    void make_invalid() {depth = MAX_RAY_REFLECTIONS;}

    Ray(vec3 i, vec3 j, int k, ivec2 l) :begin(i), dir(normalize(j)), depth(k), coords(l) {}
};



class BaseMaterial {
public:
    virtual void process(Ray &ray, vec3 pi, vec3 N, float t, vec3 col, vec3 lc, vector<vector<vec3> > &ColorMap, vector<vector<int> > &samplesCount) = 0;
    virtual ~BaseMaterial() = default;
};

class MirrorMaterial : public BaseMaterial {
public:
    void process(Ray &ray, vec3 pi, vec3 N, float t, vec3 col, vec3 lc, vector<vector<vec3> > &ColorMap, vector<vector<int> > &samplesCount) {
        if (ray.getCol() == black) {
            ray.make_invalid();
        }

        ray.reflect(pi + N * EPS, reflect(ray.getDir(), N), white);
    }
};

class DiffuseMaterial : public BaseMaterial {
public:
    void process(Ray &ray, vec3 pi, vec3 N, float t, vec3 col, vec3 lc, vector<vector<vec3> > &ColorMap, vector<vector<int> > &samplesCount) {
        if (ray.getCol() == black) {
            ray.make_invalid();
        }

        vec3 rnd = normalize(vec3(distribution1(generator), distribution1(generator), distribution1(generator)));
        if (dot(N, rnd) < 0) {
            rnd *= -1;
        }

        float dt = std::max(0.0f, dot(N, rnd));

        ray.reflect(pi + N * EPS, rnd, col * dt);
    }
};

class LightMaterial : public BaseMaterial {
public:
    void process(Ray &ray, vec3 pi, vec3 N, float t, vec3 col, vec3 lc, vector<vector<vec3> > &ColorMap, vector<vector<int> > &samplesCount) {
        ray.make_invalid();
        ray.reflect(pi, pi, col);
        ColorMap[ray.getCoords().x][ray.getCoords().y] += ray.getCol();
        ++samplesCount[ray.getCoords().x][ray.getCoords().y];
    }
};

bool planeintersect(Ray &ray, float &t, vec3 N, float D) {
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

struct Triangle {
    vec3 v0;
    vec3 v1;
    vec3 v2;
    vec3 col;
    vec3 N;
    float D;
    unique_ptr<BaseMaterial> material;
    Triangle(vec3 a, vec3 b, vec3 c, vec3 cl, unique_ptr<BaseMaterial> &&m) {
        v0 = a;
        v1 = c;
        v2 = b;
        col = cl;
        material=move(m);
        vec3 e1 = v1 - v0;
        vec3 e2 = v2 - v0;
        N = normalize(cross(e1, e2));
        D = N.x * v0.x + N.y * v0.y + N.z * v0.z;
    }

    vec3 getNormal() const {
        return N;
    }

    bool intersect(Ray &ray, float &t) const {
        float newT = INFINITY;
        if (!planeintersect(ray, newT, N, D)) {
            return false;
        }

        vec3 pi = ray.getBegin() + ray.getDir() * newT;
        float s1 = square(v0, v1, v2), s2 = square(pi, v1, v2), s3 = square(v0, pi, v2), s4 = square(v0, v1, pi);
        if (abs(s1 - s2 - s3 - s4) > EPS) {
            return false;
        }
        if (ray.getCoords().y < 100) {
            cout << 1 << endl;
        }
        t = newT;
        return true;
    }
};

void traceRay(Ray &ray, vector<vector<vec3> > &ColorMap, vector<vector<int> > &samplesCount) {

    float cube_a = 4;
    const int triangles_count = 3;
    static const Triangle triangles[triangles_count] = {Triangle(vec3(-cube_a, cube_a, 0), vec3(cube_a, cube_a, 0),
                                                        vec3(cube_a, -cube_a, 0), violet, make_unique<LightMaterial>()),

                                                        Triangle(vec3(-cube_a, cube_a, 0), vec3(cube_a, -cube_a, 0),
                                                        vec3(-cube_a, -cube_a, 0), blue, make_unique<LightMaterial>()),

                                                        Triangle(vec3(-cube_a, cube_a, -cube_a), vec3(-cube_a, cube_a, cube_a),
                                                        vec3(-cube_a, -cube_a, cube_a), red, make_unique<LightMaterial>())};


    vec3 result = black;

    float t = INFINITY;
    int i = 0, cur=-1;
    for (i = 0; i < triangles_count; ++i) {
        if (triangles[i].intersect(ray, t)) {
            cur = i;
        }
    }
    if (cur > -1) {
        vec3 pi = ray.getBegin() + ray.getDir() * t;
        vec3 N = triangles[cur].getNormal();
        triangles[cur].material->process(ray, pi, N, t, triangles[cur].col, triangles[3].col, ColorMap, samplesCount);
    } else {
        ray.make_invalid();
    }
}



int main() {
    unsigned int start_time = clock();

    vector<vector<vec3> > ColorMap(W, vector<vec3>(H, vec3(0)));
    vector<vector<int> > samplesCount(W, vector<int>(H, 0));

    bitmap_image image(H, W);

    image.clear();


    for (int i = 0; i < RAYS_PER_PIXEL; ++i) {
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                vec3 dir = vec3((x - W / 2 + distribution2(generator)) / W,
                                -(y - H / 2 + distribution2(generator)) / H,
                                1);
                Ray ray = Ray(vec3(0, 0, -20), dir, 0, ivec2(x, y));
                while (ray.is_valid()) {
                    traceRay(ray, ColorMap, samplesCount);
                }
//                if ((i + 1) % update_img == 0) {
//                    vec3 c = ColorMap[x][y] / static_cast<float>(samplesCount[x][y]) * 255.0f;
//                    image.set_pixel(x, y, c.r, c.g, c.b);
//                }
            }
        }
//        if ((i + 1) % update_img == 0) {
//            image.save_image("last.bmp");
//        }
    }

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            if (!samplesCount[x][y]) {
                continue;
            }
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