#include <memory>
#include <random>
#include <ctime>

#include "bitmap_image.hpp"
#include "glm/geometric.hpp"
#include "glm/vec3.hpp"

using namespace std;
using namespace glm;

static const vec3 white(1, 1, 1);
static const vec3 red(1, 0.1, 0.1);
static const vec3 blue(0.1, 0.1, 1);
static const vec3 violet(1, 0.1, 1);
static const vec3 yellow(1, 1, 0.1);
static const vec3 black(0, 0, 0);
static const vec3 green(0.1, 1, 0.1);

const int W = 500;
const int H = 500;
const int RAYS_PER_PIXEL = 200;
const int MAX_RAY_REFLECTIONS = 8;
const float EPS = 1e-6;

std::default_random_engine generator(time(0));
std::uniform_real_distribution<> distribution(-0.5f, 0.5f);

class Ray;

float square(vec3 A, vec3 B, vec3 C) {
    vec3 a = B - A;
    vec3 b = C - A;
    vec3 c = cross(a, b);
    return length(c) / 2;
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
    virtual void process(Ray &ray, vec3 pi, vec3 N) = 0;
    virtual ~BaseMaterial() = default;
    vec3 color;
    BaseMaterial(vec3 col) :color(col) {};
};

class MirrorMaterial : public BaseMaterial {
public:
    MirrorMaterial(vec3 col) :BaseMaterial(col) {};
    void process(Ray &ray, vec3 pi, vec3 N) {
        if (ray.getCol() == black) {
            ray.make_invalid();
            return;
        }

        ray.reflect(pi + N * EPS, reflect(ray.getDir(), N), white);
    }
};

class DiffuseMaterial : public BaseMaterial {
public:
    DiffuseMaterial(vec3 col) :BaseMaterial(col) {};
    void process(Ray &ray, vec3 pi, vec3 N) {
        if (ray.getCol() == black) {
            ray.make_invalid();
            return;
        }

        vec3 rnd = normalize(vec3(distribution(generator), distribution(generator), distribution(generator)));
        if (dot(N, rnd) < 0) {
            rnd *= -1;
        }

        float dt = std::max(0.0f, dot(N, rnd));

        ray.reflect(pi + N * EPS, rnd, color * dt);
    }
};

class LightMaterial : public BaseMaterial {
public:
    vector<vector<vec3> > &ColorMap;
    vector<vector<int> > &SamplesCount;
    LightMaterial(vector<vector<vec3> > &CM, vector<vector<int> > &SC, vec3 col) :ColorMap(CM), SamplesCount(SC), BaseMaterial(col){};
    void process(Ray &ray, vec3 pi, vec3 N) {
        ray.setCol(ray.getCol() * color);
        ColorMap[ray.getCoords().x][ray.getCoords().y] += ray.getCol();
        ++SamplesCount[ray.getCoords().x][ray.getCoords().y];
        ray.make_invalid();
    }
};

bool planeintersect(Ray &ray, float &t, vec3 N, float D) {
    vec3 o = ray.getBegin();
    vec3 d = ray.getDir();
    if (std::abs(dot(N, d)) < EPS) {
        return false;
    } else {
        float newT = -(dot(o, N) + D) / dot(d, N);
        if (t > newT && newT > 0) {
            t = newT;
            return true;
        } else {
            return false;
        }
    }
}

struct Triangle {
    std::array<vec3, 3> vertices;
    vec3 N;
    float D;
    unique_ptr<BaseMaterial> material;
    Triangle(std::array<vec3, 3> a, unique_ptr<BaseMaterial> &&m): vertices(a), material(move(m)) {
        vec3 e1 = vertices[1] - vertices[0];
        vec3 e2 = vertices[2] - vertices[0];
        N = normalize(cross(e1, e2));
        D = -dot(N, vertices[0]);
    }

    vec3 getNormal() const {
        return N;
    }

    bool intersect(Ray &ray, float &t) const {
        float newT = INFINITY;
        if (!planeintersect(ray, newT, N, D)) {
            return false;
        }

        if (newT >= t) {
            return false;
        }

        vec3 pi = ray.getBegin() + ray.getDir() * newT;
        float full_square = square(vertices[0], vertices[1], vertices[2]);
        float small_square_1 = square(pi, vertices[1], vertices[2]);
        float small_square_2 = square(vertices[0], pi, vertices[2]);
        float small_square_3 = square(vertices[0], vertices[1], pi);
        if (abs(full_square - small_square_1 - small_square_2 - small_square_3) > EPS) {
            return false;
        }
        if (newT < t) {
            t = newT;
        }
        return true;
    }
};

void traceRay(Ray &ray, const std::array<Triangle, 12> &triangles) {
    vec3 result = black;
    float t = INFINITY;
    int triangles_count = triangles.size();
    int i = 0, cur=-1;
    for (i = 0; i < triangles_count; ++i) {
        if (triangles[i].intersect(ray, t)) {
            cur = i;
        }
    }
    if (cur > -1) {
        vec3 pi = ray.getBegin() + ray.getDir() * t;
        vec3 N = triangles[cur].getNormal();
        triangles[cur].material->process(ray, pi, N);
    } else {
        ray.make_invalid();
    }
}



int main() {
    unsigned int start_time = clock();

    vector<vector<vec3> > ColorMap(W, vector<vec3>(H, vec3(0)));
    vector<vector<int> > SamplesCount(W, vector<int>(H, 0));


    float cube_a = 10;
    const int triangles_count = 12;
    const array<Triangle, triangles_count> triangles = {Triangle({vec3(-cube_a, cube_a, cube_a), vec3(cube_a, cube_a, cube_a),
                                                                 vec3(cube_a, -cube_a, cube_a)}, make_unique<DiffuseMaterial>(white)),

                                                        Triangle({vec3(-cube_a, cube_a, cube_a), vec3(cube_a, -cube_a, cube_a),
                                                                 vec3(-cube_a, -cube_a, cube_a)}, make_unique<DiffuseMaterial>(white)),

                                                        Triangle({vec3(-cube_a, cube_a, -cube_a), vec3(-cube_a, cube_a, cube_a),
                                                                 vec3(-cube_a, -cube_a, cube_a)}, make_unique<DiffuseMaterial>(red)),

                                                        Triangle({vec3(-cube_a, cube_a, -cube_a), vec3(-cube_a, -cube_a, cube_a),
                                                                 vec3(-cube_a, -cube_a, -cube_a)}, make_unique<DiffuseMaterial>(red)),

                                                        Triangle({vec3(cube_a, cube_a, -cube_a), vec3(cube_a, -cube_a, cube_a),
                                                                 vec3(cube_a, cube_a, cube_a)}, make_unique<DiffuseMaterial>(green)),

                                                        Triangle({vec3(cube_a, cube_a, -cube_a), vec3(cube_a, -cube_a, -cube_a),
                                                                 vec3(cube_a, -cube_a, cube_a)}, make_unique<DiffuseMaterial>(green)),

                                                        Triangle({vec3(-cube_a, cube_a, cube_a), vec3(cube_a, cube_a, -cube_a),
                                                                 vec3(cube_a, cube_a, cube_a)}, make_unique<DiffuseMaterial>(white)),

                                                        Triangle({vec3(-cube_a, cube_a, cube_a), vec3(-cube_a, cube_a, -cube_a),
                                                                 vec3(cube_a, cube_a, -cube_a)}, make_unique<DiffuseMaterial>(white)),

                                                        Triangle({vec3(-cube_a, -cube_a, cube_a), vec3(cube_a, -cube_a, cube_a),
                                                                 vec3(cube_a, -cube_a, -cube_a)}, make_unique<DiffuseMaterial>(white)),

                                                        Triangle({vec3(-cube_a, -cube_a, cube_a), vec3(cube_a, -cube_a, -cube_a),
                                                                 vec3(-cube_a, -cube_a, -cube_a)}, make_unique<DiffuseMaterial>(white)),

                                                        Triangle({vec3(-2, cube_a - 1, 2), vec3(2, cube_a - 1, 2),
                                                                 vec3(2, cube_a - 1, -2)}, make_unique<LightMaterial>(ColorMap, SamplesCount, white)),

                                                        Triangle({vec3(-2, cube_a - 1, 2), vec3(2, cube_a - 1, -2),
                                                                 vec3(-2, cube_a - 1, -2)}, make_unique<LightMaterial>(ColorMap, SamplesCount, white))};



    bitmap_image image(H, W);

    image.clear();


    for (int i = 0; i < RAYS_PER_PIXEL; ++i) {
        for (int y = 0; y < H; ++y) {
            #pragma omp parallel for num_threads(4)
            for (int x = 0; x < W; ++x) {
                vec3 dir = vec3((x - W / 2 + distribution(generator)) / W,
                                -(y - H / 2 + distribution(generator)) / H,
                                1);
                Ray ray = Ray(vec3(0, 0, -20), dir, 0, ivec2(x, y));
                while (ray.is_valid()) {
                    traceRay(ray, triangles);
                }
            }
        }
    }

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            if (!SamplesCount[x][y]) {
                continue;
            }
            vec3 c = ColorMap[x][y] / static_cast<float>(SamplesCount[x][y]) * 255.0f;
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