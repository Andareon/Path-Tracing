#include <memory>
#include <random>
#include <ctime>
#include <cmath>
#include <chrono>
#include <stdio.h>

#include "bitmap_image.hpp"
#include "glm/geometric.hpp"
#include "config.h"
#include "ray.h"

using namespace std;
using namespace glm;

const float PI = 3.141593;

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
    vector<vector<vec3> > Ans(Config::get().width, vector<vec3>(Config::get().height, vec3(0, 0, 0)));
    for (int i = 0; i < Config::get().height; ++i) {
        for (int j = 0; j < Config::get().width; ++j) {
            vec3 val = vec3(0, 0, 0);
            float wsum = 0;
            for (int iy = i - rs; iy <= i + rs; ++iy) {
                for (int ix = j - rs; ix <= j + rs; ++ix) {
                    int x = std::min(Config::get().width - 1, std::max(0, ix));
                    int y = std::min(Config::get().height- 1, std::max(0, iy));
                    int dsq = (ix - j) * (ix - j) + (iy - i) * (iy - i);
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

vector<vector<vec3> > median_filter(vector<vector<vec3> > ColorMap, int t) {
    vector<vector<vec3> > Ans(Config::get().width, vector<vec3>(Config::get().height, vec3(0, 0, 0)));
    for (int y = 0; y < Config::get().height; ++y) {
        for (int x = 0; x < Config::get().width; ++x) {
            vector<float> win_r;
            vector<float> win_g;
            vector<float> win_b;
            for (int ii = -t; ii < t + 1; ++ii) {
                for (int jj = -t; jj < t + 1; ++jj) {
                    int i = std::max(std::min(ii + x, Config::get().width - 1), 0);
                    int j = std::max(std::min(jj + y, Config::get().height - 1), 0);
                    win_r.push_back(ColorMap[i][j].r);
                    win_g.push_back(ColorMap[i][j].g);
                    win_b.push_back(ColorMap[i][j].b);
                }
            }
            sort(win_r.begin(), win_r.end());
            Ans[x][y].r = win_r[t * t / 2];

            sort(win_g.begin(), win_g.end());
            Ans[x][y].g = win_g[t * t / 2];

            sort(win_b.begin(), win_b.end());
            Ans[x][y].b = win_b[t * t / 2];
        }
    }
    return Ans;
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
        ray.reflect(pi + N * Config::get().EPS, reflect(ray.getDir(), N), color);
    }
};

class DiffuseMaterial : public BaseMaterial {
public:
    DiffuseMaterial(vec3 col) :BaseMaterial(col) {};
    void process(Ray &ray, vec4 pi, vec4 N) {
        static default_random_engine generator(time(0));
        static uniform_real_distribution<> distribution(.0f, 1.f);
        float xi1 = distribution(generator),
              xi2 = distribution(generator);
        vec4 rnd = normalize(vec4(sqrt(xi1) * cos(2 * PI * xi2), sqrt(xi1) * sin(2 * PI * xi2), sqrt(1 - xi1), 0));
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
    vector<vector<vec3> > &Color2Map;
    vector<vector<int> > &SamplesCount;
public:
    LightMaterial(vector<vector<vec3> > &CM, vector<vector<vec3> > &C2M, vector<vector<int> > &SC, vec3 col) :ColorMap(CM), Color2Map(C2M), SamplesCount(SC), BaseMaterial(col){};
    void process(Ray &ray, vec4 pi, vec4 N) {
        if (dot(ray.getDir(), N) < 0) {
            ray.make_invalid();
            return;
        }
        ColorMap[ray.getCoords().x][ray.getCoords().y] += ray.getCol() * color;
        Color2Map[ray.getCoords().x][ray.getCoords().y] += (ray.getCol() * color) * (ray.getCol() * color);
        ++SamplesCount[ray.getCoords().x][ray.getCoords().y];
        ray.make_invalid();
    }
};

class TransparentMaterial : public BaseMaterial {
private:
    float ref_in;
public:
    TransparentMaterial(vec3 col, float n) :BaseMaterial(col), ref_in(n){};
    void process(Ray &ray, vec4 pi, vec4 N) {
        vec4 dir;
        vec3 col;
        if (dot(N, ray.getDir()) > 0) {
            dir = refract(ray.getDir(), -N, ref_in);
            col = vec3(1, 1, 1);
        } else {
            dir = refract(ray.getDir(), N, 1 / ref_in);
            col = color;
        }

        ray.reflect(pi + dir * Config::get().EPS, dir, col);
    }
};

bool planeintersect(Ray &ray, float &t, vec4 plane) {
    vec4 o = ray.getBegin();
    vec4 d = ray.getDir();
    if (abs(dot(plane, d)) < Config::get().EPS) {
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
    vec4 plane;
    BaseMaterial* material;
    array<vec4, 3> vertices;

public:
    Triangle(array<vec4, 3> a, BaseMaterial* m): vertices(a), material(move(m)) {
        vec4 e1 = vertices[1] - vertices[0];
        vec4 e2 = vertices[2] - vertices[0];
        plane = normalize(cross(e1, e2));
        plane.w = -dot(plane, vertices[0]);
    }

    vec4 getNormal() const {
        return vec4(plane.x, plane.y, plane.z, 0);
    }

    void setNormal(vec3 N) {
        N = normalize(N);
        plane.x = N.x;
        plane.y = N.y;
        plane.z = N.z;
        plane.w = 0;
        plane.w = -dot(plane, vertices[0]);
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

class Scene {
private:
    vector<Triangle> triangles;
public:
    void LoadModel(const char *path, BaseMaterial *material) {
        vector<vec4> temp_vertices;
        vector<vec2> temp_uvs;
        vector<vec3> temp_normals;
        FILE *file = fopen(path, "r");
        if (file == NULL) {
            printf("Impossible to open the file !\n");
            return;
        }
        while (1) {
            char lineHeader[128];
            int res = fscanf(file, "%s", lineHeader);
            if (res == EOF)
                break;
            if (strcmp(lineHeader, "v") == 0) {
                glm::vec4 vertex = vec4(1);
                fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z);
                temp_vertices.push_back(vertex);
            } else if (strcmp(lineHeader, "vt") == 0) {
                glm::vec2 uv;
                fscanf(file, "%f %f\n", &uv.x, &uv.y);
                temp_uvs.push_back(uv);
            } else if (strcmp(lineHeader, "vn") == 0) {
                glm::vec3 normal;
                fscanf(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z);
                temp_normals.push_back(normal);
            } else if (strcmp(lineHeader, "f") == 0) {
                int vi1, vi2, vi3, uvi1, uvi2, uvi3, ni1, ni2, ni3;
                int matches = fscanf(file, " %d/%d/%d %d/%d/%d %d/%d/%d\n", &vi1, &uvi1, &ni1, &vi2, &uvi2, &ni2, &vi3,
                                     &uvi3, &ni3);
                vi1--;
                vi2--;
                vi3--;
                triangles.push_back(Triangle({temp_vertices[vi1], temp_vertices[vi2], temp_vertices[vi3]}, material));
                triangles.back().setNormal(temp_normals[ni1 - 1]);
            }
        }
    }

    void AddTriangle(Triangle tr) {
        triangles.push_back(tr);
    }


    void traceRay(Ray &ray) {
        float t = INFINITY;
        int triangles_count = triangles.size();
        int i = 0, cur = -1;
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
};

int main(int argc, char* argv[]) {

    chrono::milliseconds start_time = chrono::duration_cast<chrono::milliseconds>(
            chrono::system_clock::now().time_since_epoch());

    Config::get().set_config(argc, argv);
    static default_random_engine generator(time(0));
    static uniform_real_distribution<> distribution(-0.5f, 0.5f);


    vector<vector<vec3> > ColorMap(Config::get().width, vector<vec3>(Config::get().height, vec3(0, 0, 0)));
    vector<vector<vec3> > Color2Map(Config::get().width, vector<vec3>(Config::get().height, vec3(0, 0, 0)));
    vector<vector<int> > SamplesCount(Config::get().width, vector<int>(Config::get().height, 0));
    vector<BaseMaterial*> Materials = {new DiffuseMaterial(vec3(1, 1, 1)),
                                       new DiffuseMaterial(vec3(1, 0, 0)),
                                       new DiffuseMaterial(vec3(0, 1, 0)),
                                       new LightMaterial(ColorMap, Color2Map, SamplesCount, vec3(1, 1, 1)),
                                       new TransparentMaterial(vec3(1, 1, 0), 1.25),
                                       new MirrorMaterial(vec3(1, 1, 1))};

    Scene scene;
    scene.LoadModel("../1.obj", Materials[2]);

    float cube_a = 10;
    scene.AddTriangle(Triangle({vec4(-cube_a, cube_a, cube_a, 1), vec4(cube_a, cube_a, cube_a, 1),
                                vec4(cube_a, -cube_a, cube_a, 1)}, Materials[0]));

    scene.AddTriangle(Triangle({vec4(-cube_a, cube_a, cube_a, 1), vec4(cube_a, -cube_a, cube_a, 1),
                                vec4(-cube_a, -cube_a, cube_a, 1)}, Materials[0]));

    scene.AddTriangle(Triangle({vec4(-cube_a, cube_a, -cube_a, 1), vec4(-cube_a, cube_a, cube_a, 1),
                                vec4(-cube_a, -cube_a, cube_a, 1)}, Materials[1]));

    scene.AddTriangle(Triangle({vec4(-cube_a, cube_a, -cube_a, 1), vec4(-cube_a, -cube_a, cube_a, 1),
                                vec4(-cube_a, -cube_a, -cube_a, 1)}, Materials[1]));

    scene.AddTriangle(Triangle({vec4(cube_a, cube_a, -cube_a, 1), vec4(cube_a, -cube_a, cube_a, 1),
                                vec4(cube_a, cube_a, cube_a, 1)}, Materials[2]));

    scene.AddTriangle(Triangle({vec4(cube_a, cube_a, -cube_a, 1), vec4(cube_a, -cube_a, -cube_a, 1),
                                vec4(cube_a, -cube_a, cube_a, 1)}, Materials[2]));

    scene.AddTriangle(Triangle({vec4(-cube_a, cube_a, cube_a, 1), vec4(cube_a, cube_a, -cube_a, 1),
                                vec4(cube_a, cube_a, cube_a, 1)}, Materials[0]));

    scene.AddTriangle(Triangle({vec4(-cube_a, cube_a, cube_a, 1), vec4(-cube_a, cube_a, -cube_a, 1),
                                vec4(cube_a, cube_a, -cube_a, 1)}, Materials[0]));

    scene.AddTriangle(Triangle({vec4(-cube_a, -cube_a, cube_a, 1), vec4(cube_a, -cube_a, cube_a, 1),
                                vec4(cube_a, -cube_a, -cube_a, 1)}, Materials[0]));

    scene.AddTriangle(Triangle({vec4(-cube_a, -cube_a, cube_a, 1), vec4(cube_a, -cube_a, -cube_a, 1),
                                vec4(-cube_a, -cube_a, -cube_a, 1)}, Materials[0]));

    scene.AddTriangle(Triangle({vec4(-1, cube_a - 1, 1, 1), vec4(1, cube_a - 1, 1, 1),
                                vec4(1, cube_a - 1, -1, 1)}, Materials[3]));

    scene.AddTriangle(Triangle({vec4(-1, cube_a - 1, 1, 1), vec4(1, cube_a - 1, -1, 1),
                                vec4(-1, cube_a - 1, -1, 1)}, Materials[3]));




    bitmap_image image(Config::get().width, Config::get().height);

    image.clear();

    for (int i = 0; i < Config::get().RAYS_PER_PIXEL; ++i) {
        for (int y = 0; y < Config::get().height; ++y) {
            #pragma omp parallel for num_threads(4)
            for (int x = 0; x < Config::get().width; ++x) {
                float SampleCount = static_cast<float>(SamplesCount[x][y]);
                if (i > 10 && SampleCount) {
                    vec3 D = (Color2Map[x][y] / SampleCount - (ColorMap[x][y] / SampleCount) * (ColorMap[x][y] / SampleCount));
                    if ((i % 4) && D.r < Config::get().error && D.g < Config::get().error && D.b < Config::get().error) {
                        continue;
                    }
                }
                vec4 dir = vec4((static_cast<float>(x) - Config::get().width / 2 + distribution(generator)) / Config::get().width,
                                -(static_cast<float>(y) - Config::get().height / 2 + distribution(generator)) / Config::get().height,
                                1, 0);
                Ray ray = Ray(vec4(0, 0, -20, 1), dir, 0, ivec2(x, y));
                while (ray.is_valid()) {
                    scene.traceRay(ray);
                }
                if (i % Config::get().update == 0) {
                    if (SamplesCount[x][y]) {
                        vec3 c = pow(ColorMap[x][y] / SampleCount, vec3(1/2.2)) * 255.0f;
                        image.set_pixel(x, y, c.r, c.g, c.b);
                    }
                }
            }
        }
        if (i % Config::get().update == 0) {
            image.save_image("result.bmp");
            cerr << "Image update" << endl;
        }
        cerr << i + 1 << " rays per pixel were sent" << endl;
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
    if (Config::get().median) {
        ColorMap = median_filter(ColorMap, Config::get().median);
    }
    for (int y = 0; y < Config::get().height; ++y) {
        for (int x = 0; x < Config::get().width; ++x) {
            if (!SamplesCount[x][y]) {
                continue;
            }
            vec3 c = ColorMap[x][y];
            image.set_pixel(x, y, c.r, c.g, c.b);
        }
    }

    chrono::milliseconds end_time = chrono::duration_cast<chrono::milliseconds>(
            chrono::system_clock::now().time_since_epoch());
    time_t t = time(0);
    struct tm * now = localtime( & t );
    string date = to_string(now->tm_year + 1900) + '-' + to_string(now->tm_mon + 1) + '-' +
                  to_string(now->tm_mday) + '-' + to_string(now->tm_hour) + '-' + to_string(now->tm_min) + '-' +
                  to_string(now->tm_sec) + "  " + to_string((end_time - start_time).count()) + "   " + to_string(Config::get().RAYS_PER_PIXEL);
    image.save_image(date + ".bmp");
    image.save_image("result.bmp");
}