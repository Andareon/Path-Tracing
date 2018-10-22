#include <memory>
#include <random>
#include <ctime>
#include <cmath>

#include "bitmap_image.hpp"
#include "glm/geometric.hpp"
#include "config.h"
#include "ray.h"


using namespace std;
using namespace glm;

static const vec3 white(1, 1, 1);
static const vec3 red(1, 0.1, 0.1);
static const vec3 blue(0.1, 0.1, 1);
static const vec3 violet(1, 0.1, 1);
static const vec3 yellow(1, 1, 0.1);
static const vec3 black(0, 0, 0);
static const vec3 green(0.1, 1, 0.1);
static const vec3 gray(0.5, 0.5, 0.5);

const float PI = 3.141593;

std::default_random_engine generator(time(0));
std::uniform_real_distribution<> distribution(-0.5f, 0.5f);

float square(vec4 A, vec4 B, vec4 C) {
    vec3 a = B - A;
    vec3 b = C - A;
    vec3 c = cross(a, b);
    return length(c) / 2;
}

vec4 cross(vec4 a, vec4 b) {
    return vec4(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x, 0);
}

vector<vector<vec3> > gauss_blur(vector<vector<vec3> > ColorMap, float r) {
    int rs = ceil(r * 2.57);
    vector<vector<vec3> > Ans(Config::get().width, vector<vec3>(Config::get().height, vec3(0)));
    for (int i = 0; i < Config::get().height; ++i) {
        for (int j = 0; j < Config::get().width; ++j) {
            vec3 val = vec3(0);
            float wsum = 0;
            for (int iy = i - rs; iy <= i + rs; ++iy) {
                for (int ix = j - rs; ix <= j + rs; ++ix) {
                    int x = std::min(Config::get().width-1, std::max(0, ix));
                    int y = std::min(Config::get().height-1, std::max(0, iy));
                    int dsq = (ix-j)*(ix-j)+(iy-i)*(iy-i);
                    float wght = exp(-dsq / (2 * r * r)) / (PI * 2 * r * r);
                    val += ColorMap[x][y] * wght;
                    wsum += wght;
                }
            }
            Ans[j][i] = round(val / wsum);
        }
    }
    return Ans;
}

vec4 refract(vec4 dir, vec4 N, float n) {
    float cosI = -dot(N, dir);
    float sinT = n * sqrt(1.0f - cosI * cosI);
    if (sinT > 1) {
        return reflect(dir, N);
    }
    float cosT = sqrt(1.0f - sinT * sinT);
    vec4 dirT = dir * n + N * (n * cosI - cosT);
    return dirT;
}



class BaseMaterial {
protected:
    vec3 color;
public:
    virtual void process(Ray &ray, vec4 pi, vec4 N) = 0;
    virtual ~BaseMaterial() = default;
    BaseMaterial(vec3 col) :color(col) {};
};

class MirrorMaterial : public BaseMaterial {
public:
    MirrorMaterial(vec3 col) :BaseMaterial(col) {};
    void process(Ray &ray, vec4 pi, vec4 N) {
        if (ray.getCol() == black) {
            ray.make_invalid();
            return;
        }

        ray.reflect(pi + N * Config::get().EPS, reflect(ray.getDir(), N), color);
    }
};

class DiffuseMaterial : public BaseMaterial {
public:
    DiffuseMaterial(vec3 col) :BaseMaterial(col) {};
    void process(Ray &ray, vec4 pi, vec4 N) {
        if (ray.getCol() == black) {
            ray.make_invalid();
            return;
        }

        vec4 rnd = normalize(vec4(distribution(generator), distribution(generator), distribution(generator), 0));
        if (dot(N, rnd) < 0) {
            rnd *= -1;
        }

        float dt = std::max(0.0f, dot(N, rnd));

        ray.reflect(pi + N * Config::get().EPS, rnd, color * dt);
    }
};

class LightMaterial : public BaseMaterial {
private:
    vector<vector<vec3> > &ColorMap;
    vector<vector<int> > &SamplesCount;
public:
    LightMaterial(vector<vector<vec3> > &CM, vector<vector<int> > &SC, vec3 col) :ColorMap(CM), SamplesCount(SC), BaseMaterial(col){};
    void process(Ray &ray, vec4 pi, vec4 N) {
        if (dot(ray.getDir(), N) < 0) {
            ray.make_invalid();
            return;
        }
        ColorMap[ray.getCoords().x][ray.getCoords().y] += ray.getCol() * color;
        ++SamplesCount[ray.getCoords().x][ray.getCoords().y];
        ray.make_invalid();
    }
};

class TransparentMaterial : public BaseMaterial {
private:
    float n1;
public:
    TransparentMaterial(vec3 col, float n) :BaseMaterial(col), n1(n){};
    void process(Ray &ray, vec4 pi, vec4 N) {
        if (ray.getCol() == black) {
            ray.make_invalid();
            return;
        }

        vec4 dir;
        vec3 col;
        if (dot(N, ray.getDir()) > 0) {
            dir = refract(ray.getDir(), -N, n1);
            col = white;
        } else {
            dir = refract(ray.getDir(), N, 1 / n1);
            col = color;
        }

        ray.reflect(pi + dir * Config::get().EPS, dir, col);
    }
};


bool planeintersect(Ray &ray, float &t, vec4 plane) {
    vec4 o = ray.getBegin();
    vec4 d = ray.getDir();
    if (std::abs(dot(plane, d)) < Config::get().EPS) {
        return false;
    } else {
        float newT = -dot(o, plane) / dot(d, plane);
        if (t > newT && newT > 0) {
            t = newT;
            return true;
        } else {
            return false;
        }
    }
}

class Triangle {
private:
    std::array<vec4, 3> vertices;
    vec4 plane;
    unique_ptr<BaseMaterial> material;

public:
    Triangle(std::array<vec4, 3> a, unique_ptr<BaseMaterial> &&m): vertices(a), material(move(m)) {
        vec4 e1 = vertices[1] - vertices[0];
        vec4 e2 = vertices[2] - vertices[0];
        plane = normalize(cross(e1, e2));
        plane.w = -dot(plane, vertices[0]);
    }

    vec4 getNormal() const {
        return vec4(plane.x, plane.y, plane.z, 0);
    }

    BaseMaterial &getMaterial() const {
        return *material;
    }

    bool intersect(Ray &ray, float &t) const {
        float newT = INFINITY;
        if (!planeintersect(ray, newT, plane)) {
            return false;
        }

        if (newT >= t) {
            return false;
        }

        vec4 pi = ray.getBegin() + ray.getDir() * newT;
        float full_square = square(vertices[0], vertices[1], vertices[2]);
        float small_square_1 = square(pi, vertices[1], vertices[2]);
        float small_square_2 = square(vertices[0], pi, vertices[2]);
        float small_square_3 = square(vertices[0], vertices[1], pi);
        if (abs(full_square - small_square_1 - small_square_2 - small_square_3) > Config::get().EPS) {
            return false;
        }
            t = newT;
        return true;
    }
};

void traceRay(Ray &ray, const std::array<Triangle, 18> &triangles) {
    float t = INFINITY;
    int triangles_count = triangles.size();
    int i = 0, cur=-1;
    for (i = 0; i < triangles_count; ++i) {
        if (triangles[i].intersect(ray, t)) {
            cur = i;
        }
    }
    if (cur > -1) {
        vec4 pi = ray.getBegin() + ray.getDir() * t;
        vec4 N = triangles[cur].getNormal();
        triangles[cur].getMaterial().process(ray, pi, N);
    } else {
        ray.make_invalid();
    }
}



int main(int argc, char* argv[]) {
    unsigned int start_time = clock();

    Config::get().set_config(argc, argv);

    vector<vector<vec3> > ColorMap(Config::get().width, vector<vec3>(Config::get().height, vec3(0)));
    vector<vector<int> > SamplesCount(Config::get().width, vector<int>(Config::get().height, 0));

    float cube_a = 10;
    float pir_a = 4;
    const int triangles_count = 18;
    const array<Triangle, triangles_count> triangles = {Triangle({vec4(-cube_a, cube_a, cube_a, 1), vec4(cube_a, cube_a, cube_a, 1),
                                                                 vec4(cube_a, -cube_a, cube_a, 1)}, make_unique<DiffuseMaterial>(white)),

                                                        Triangle({vec4(-cube_a, cube_a, cube_a, 1), vec4(cube_a, -cube_a, cube_a, 1),
                                                                 vec4(-cube_a, -cube_a, cube_a, 1)}, make_unique<DiffuseMaterial>(white)),

                                                        Triangle({vec4(-cube_a, cube_a, -cube_a, 1), vec4(-cube_a, cube_a, cube_a, 1),
                                                                 vec4(-cube_a, -cube_a, cube_a, 1)}, make_unique<DiffuseMaterial>(red)),

                                                        Triangle({vec4(-cube_a, cube_a, -cube_a, 1), vec4(-cube_a, -cube_a, cube_a, 1),
                                                                 vec4(-cube_a, -cube_a, -cube_a, 1)}, make_unique<DiffuseMaterial>(red)),

                                                        Triangle({vec4(cube_a, cube_a, -cube_a, 1), vec4(cube_a, -cube_a, cube_a, 1),
                                                                 vec4(cube_a, cube_a, cube_a, 1)}, make_unique<DiffuseMaterial>(green)),

                                                        Triangle({vec4(cube_a, cube_a, -cube_a, 1), vec4(cube_a, -cube_a, -cube_a, 1),
                                                                 vec4(cube_a, -cube_a, cube_a, 1)}, make_unique<DiffuseMaterial>(green)),

                                                        Triangle({vec4(-cube_a, cube_a, cube_a, 1), vec4(cube_a, cube_a, -cube_a, 1),
                                                                 vec4(cube_a, cube_a, cube_a, 1)}, make_unique<DiffuseMaterial>(white)),

                                                        Triangle({vec4(-cube_a, cube_a, cube_a, 1), vec4(-cube_a, cube_a, -cube_a, 1),
                                                                 vec4(cube_a, cube_a, -cube_a, 1)}, make_unique<DiffuseMaterial>(white)),

                                                        Triangle({vec4(-cube_a, -cube_a, cube_a, 1), vec4(cube_a, -cube_a, cube_a, 1),
                                                                 vec4(cube_a, -cube_a, -cube_a, 1)}, make_unique<DiffuseMaterial>(white)),

                                                        Triangle({vec4(-cube_a, -cube_a, cube_a, 1), vec4(cube_a, -cube_a, -cube_a, 1),
                                                                 vec4(-cube_a, -cube_a, -cube_a, 1)}, make_unique<DiffuseMaterial>(white)),

                                                        Triangle({vec4(-2, cube_a - 1, 2, 1), vec4(2, cube_a - 1, 2, 1),
                                                                 vec4(2, cube_a - 1, -2, 1)}, make_unique<LightMaterial>(ColorMap, SamplesCount, white)),

                                                        Triangle({vec4(-2, cube_a - 1, 2, 1), vec4(2, cube_a - 1, -2, 1),
                                                                 vec4(-2, cube_a - 1, -2, 1)}, make_unique<LightMaterial>(ColorMap, SamplesCount, white)),

                                                        Triangle({vec4(5 + 0, pir_a, 0, 1), vec4(5 + -pir_a, -1, pir_a, 1),
                                                                  vec4(5 + pir_a, -1, pir_a, 1)}, make_unique<TransparentMaterial>(yellow, 1.25)),

                                                        Triangle({vec4(5 + 0, pir_a, 0, 1), vec4(5 + pir_a, -1, pir_a, 1),
                                                                  vec4(5 + pir_a, -1, -pir_a, 1)}, make_unique<TransparentMaterial>(yellow, 1.25)),

                                                        Triangle({vec4(5 + 0, pir_a, 0, 1), vec4(5 + -pir_a, -1, -pir_a, 1),
                                                                  vec4(5 + -pir_a, -1, pir_a, 1)}, make_unique<TransparentMaterial>(yellow, 1.25)),

                                                        Triangle({vec4(5 + 0, pir_a, 0, 1), vec4(5 + pir_a, -1, -pir_a, 1),
                                                                  vec4(5 + -pir_a, -1, -pir_a, 1)}, make_unique<TransparentMaterial>(yellow, 1.25)),

                                                        Triangle({vec4(5 + -pir_a, -1, pir_a, 1), vec4(5 + pir_a, -1, -pir_a, 1),
                                                                  vec4(5 + pir_a, -1, pir_a, 1)}, make_unique<TransparentMaterial>(yellow, 1.25)),

                                                        Triangle({vec4(5 + -pir_a, -1, pir_a, 1), vec4(5 + -pir_a, -1, -pir_a, 1),
                                                                  vec4(5 + pir_a, -1, -pir_a, 1)}, make_unique<TransparentMaterial>(yellow, 1.25))};



    bitmap_image image(Config::get().width, Config::get().height);

    image.clear();


    for (int i = 0; i < Config::get().RAYS_PER_PIXEL; ++i) {
        for (int y = 0; y < Config::get().height; ++y) {
            #pragma omp parallel for num_threads(4)
            for (int x = 0; x < Config::get().width; ++x) {
                vec4 dir = vec4((x - Config::get().width / 2 + distribution(generator)) / Config::get().width,
                                -(y - Config::get().height / 2 + distribution(generator)) / Config::get().height,
                                1, 0);
                Ray ray = Ray(vec4(0, 0, -20, 1), dir, 0, ivec2(x, y));
                while (ray.is_valid()) {
                    traceRay(ray, triangles);
                }
                if ((i + 1) % 30 == 0) {
                    if (SamplesCount[x][y]) {
                        vec3 c = pow(ColorMap[x][y] / static_cast<float>(SamplesCount[x][y]), vec3(1/2.2)) * 255.0f;
                        image.set_pixel(x, y, c.r, c.g, c.b);
                    }
                }
            }
        }
        if ((i + 1) % 30 == 0) {
            image.save_image("result.bmp");
        }
    }

    for (int y = 0; y < Config::get().height; ++y) {
        for (int x = 0; x < Config::get().width; ++x) {
            if (!SamplesCount[x][y]) {
                continue;
            }
            ColorMap[x][y] = pow(ColorMap[x][y] / static_cast<float>(SamplesCount[x][y]), vec3(1/2.2)) * 255.0f;
        }
    }

//    ColorMap = gauss_blur(ColorMap, 0);
    for (int y = 0; y < Config::get().height; ++y) {
        for (int x = 0; x < Config::get().width; ++x) {
            if (!SamplesCount[x][y]) {
                continue;
            }
            vec3 c = ColorMap[x][y];
            image.set_pixel(x, y, c.r, c.g, c.b);
        }
    }

    unsigned int end_time = clock();
    time_t t = time(0);
    struct tm * now = localtime( & t );
    string date = to_string(now->tm_year + 1900) + '-' + to_string(now->tm_mon + 1) + '-' +
                  to_string(now->tm_mday) + '-' + to_string(now->tm_hour) + '-' + to_string(now->tm_min) + '-' +
                  to_string(now->tm_sec) + "  " + to_string(end_time - start_time) + "   " + to_string(Config::get().RAYS_PER_PIXEL);
    image.save_image(date + ".bmp");
    image.save_image("result.bmp");
}