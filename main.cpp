#include <chrono>

#include "bitmap_image.hpp"

#include "material.h"

using namespace std;
using namespace glm;

bool PlaneIntersect(Ray &ray, float &distance, vec4 plane) {
    vec4 begin = ray.GetBegin();
    vec4 direction = ray.GetDirection();
    if (abs(dot(plane, direction)) < Config::get().eps) {
        return false;
    } else {
        float new_distance = -dot(begin, plane) / dot(direction, plane);
        if (distance > new_distance && new_distance > 0) {
            distance = new_distance;
            return true;
        } else {
            return false;
        }
    }
}

vec4 Cross(vec4 a, vec4 b) {
    return vec4(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x, 0);
}

float Square(vec4 A, vec4 B, vec4 C) {
    vec4 AB = B - A;
    vec4 AC = C - A;
    vec4 cross = Cross(AB, AC);
    return length(cross) / 2;
}

vector<string> Split(string str, char sep) {
    vector<string> ans;
    int begin, end;
    for (begin = 0, end = str.find(sep); end != string::npos;
         begin = end + 1, end = str.find(sep, begin)) {
        ans.push_back(str.substr(begin, end - begin));
    }
    ans.push_back(str.substr(begin));
    return ans;
}

vector<vector<vec3> > GaussBlur(vector<vector<vec3> > color_map, float r) {
    int rs = ceil(r * 2.57);
    vector<vector<vec3> > ans(Config::get().width,
                              vector<vec3>(Config::get().height, vec3(0)));
    for (int i = 0; i < Config::get().height; ++i) {
        for (int j = 0; j < Config::get().width; ++j) {
            vec3 val = vec3(0);
            float wsum = 0;
            for (int iy = i - rs; iy <= i + rs; ++iy) {
                for (int ix = j - rs; ix <= j + rs; ++ix) {
                    int x = std::min(Config::get().width - 1, std::max(0, ix));
                    int y = std::min(Config::get().height - 1, std::max(0, iy));
                    int dsq = (ix - j) * (ix - j) + (iy - i) * (iy - i);
                    float wght = exp(-dsq / (2 * r * r)) / (pi * 2 * r * r);
                    val += color_map[x][y] * wght;
                    wsum += wght;
                }
            }
            ans[j][i] = round(val / wsum);
        }
    }
    return ans;
}

vec4 Refract(vec4 incidence_direction, vec4 normal, float refractive) {
    float cos_incidence = -dot(normal, incidence_direction);
    float sin_refraction =
            refractive * sqrt(1.0f - cos_incidence * cos_incidence);
    if (sin_refraction > 1) {
        return reflect(incidence_direction, normal);
    }
    float cos_refraction = sqrt(1.0f - sin_refraction * sin_refraction);
    vec4 refraction_direction =
            incidence_direction * refractive +
            normal * (refractive * cos_incidence - cos_refraction);
    return refraction_direction;
}

vector<vector<vec3> > MedianFilter(vector<vector<vec3> > color_map,
                                   int window_size) {
    vector<vector<vec3> > ans(Config::get().width,
                              vector<vec3>(Config::get().height, vec3(0)));
    for (int y = 0; y < Config::get().height; ++y) {
        for (int x = 0; x < Config::get().width; ++x) {
            vector<float> window_r;
            vector<float> window_g;
            vector<float> window_b;
            for (int window_x = -window_size; window_x < window_size + 1;
                 ++window_x) {
                for (int window_y = -window_size; window_y < window_size + 1;
                     ++window_y) {
                    int i = std::max(std::min(window_x + x, Config::get().width - 1), 0);
                    int j = std::max(std::min(window_y + y, Config::get().height - 1), 0);
                    window_r.push_back(color_map[i][j].r);
                    window_g.push_back(color_map[i][j].g);
                    window_b.push_back(color_map[i][j].b);
                }
            }
            sort(window_r.begin(), window_r.end());
            ans[x][y].r = window_r[window_size * window_size / 2];

            sort(window_g.begin(), window_g.end());
            ans[x][y].g = window_g[window_size * window_size / 2];

            sort(window_b.begin(), window_b.end());
            ans[x][y].b = window_b[window_size * window_size / 2];
        }
    }
    return ans;
}

class Triangle {
private:
    vec4 plane_;
    Material material_;
    array<vec4, 3> vertices_;

public:
    Triangle(const array<vec4, 3> &vertices, Material material)
            : vertices_(vertices), material_(material) {
        vec4 AB = vertices_[1] - vertices_[0];
        vec4 AC = vertices_[2] - vertices_[0];
        plane_ = normalize(Cross(AB, AC));
        plane_.w = -dot(plane_, vertices_[0]);
    }

    vec4 GetNormal() const { return vec4(plane_.x, plane_.y, plane_.z, 0); }

    void SetNormal(vec3 normal) {
        normal = normalize(normal);
        plane_.x = normal.x;
        plane_.y = normal.y;
        plane_.z = normal.z;
        plane_.w = 0;
        plane_.w = -dot(plane_, vertices_[0]);
    }

    Material GetMaterial() const { return material_; }

    bool Intersect(Ray &ray, float &distance) const {
        float new_distance = INFINITY;
        if (!PlaneIntersect(ray, new_distance, plane_)) {
            return false;
        }

        if (new_distance >= distance) {
            return false;
        }

        vec4 drop_point = ray.GetBegin() + ray.GetDirection() * new_distance;
        float full_square = Square(vertices_[0], vertices_[1], vertices_[2]);
        float small_square_1 = Square(drop_point, vertices_[1], vertices_[2]);
        float small_square_2 = Square(vertices_[0], drop_point, vertices_[2]);
        float small_square_3 = Square(vertices_[0], vertices_[1], drop_point);
        if (abs(full_square - small_square_1 - small_square_2 - small_square_3) >
            Config::get().eps) {
            return false;
        }
        distance = new_distance;
        return true;
    }
};

class Scene {
private:
    vector<Triangle> triangles_;
    vector<vector<vec3> > &color_map_;
    vector<vector<vec3> > &color2_map_;
    vector<vector<int> > &samples_count_;
    vector<Material> materials_;
    bitmap_image skybox_;

public:
    Scene(vector<vector<vec3> > &color_map, vector<vector<vec3> > &color2_map,
          vector<vector<int> > &samples_count)
            : color_map_(color_map),
              color2_map_(color2_map),
              samples_count_(samples_count) {
        if (Config::get().skybox != "") {
            skybox_ = bitmap_image(Config::get().skybox);
        }
    };
    void LoadModel(string path) {
        vector<vec4> temp_vertices;
        vector<vec2> temp_texture_coords;
        vector<vec3> temp_normals;
        ifstream file(path);
        string input;
        if (file.eof()) {
            cout << "Impossible to open the file !" << endl;
            return;
        }
        int current_material = 0;
        while (!file.eof()) {
            file >> input;
            if (input == "mtllib") {
                string mtl_file_name;
                file >> mtl_file_name;
                ifstream mtl_file(Config::get().model_path + mtl_file_name);

                string mtl_input = "1";
                MaterialCharacteristics characteristics;
                while (!mtl_file.eof()) {
                    while (!mtl_file.eof() && mtl_input != "newmtl") {
                        mtl_file >> mtl_input;
                    }
                    mtl_file >> mtl_input;
                    while (!mtl_file.eof() && mtl_input != "newmtl") {
                        if (mtl_input == "Kd") {
                            mtl_file >> characteristics.Kd.x >> characteristics.Kd.y >>
                                     characteristics.Kd.z;
                        } else if (mtl_input == "Ke") {
                            mtl_file >> characteristics.Ke.x >> characteristics.Ke.y >>
                                     characteristics.Ke.z;
                        }
                        mtl_file >> mtl_input;
                    }
                    materials_.emplace_back(Factory(characteristics, color_map_,
                                                    color2_map_, samples_count_));
                }

            } else if (input == "v") {
                vec4 vertex = vec4(1);
                file >> vertex.x >> vertex.y >> vertex.z;
                temp_vertices.push_back(vertex);
            } else if (input == "vt") {
                vec2 texture_coord;
                file >> texture_coord.x >> texture_coord.y;
                temp_texture_coords.push_back(texture_coord);
            } else if (input == "vn") {
                vec3 normal;
                file >> normal.x >> normal.y >> normal.z;
                temp_normals.push_back(normal);
            } else if (input == "f") {
                int vertex_index[3];
                int texture_coords_index[3];
                int normal_index[3];
                for (int i = 0; i < 3; ++i) {
                    string cr;
                    file >> cr;
                    vector<string> cur = Split(cr, '/');
                    cur.resize(3);
                    vertex_index[i] = atoi(cur[0].c_str()) - 1;
                    texture_coords_index[i] = atoi(cur[1].c_str()) - 1;
                    normal_index[i] = atoi(cur[2].c_str()) - 1;
                }
                triangles_.emplace_back(array<vec4, 3>{temp_vertices[vertex_index[0]],
                                                       temp_vertices[vertex_index[1]],
                                                       temp_vertices[vertex_index[2]]},
                                        materials_[current_material]);
                if (normal_index[0] >= 0) {
                    triangles_.back().SetNormal(temp_normals[normal_index[0]]);
                }
            } else if (input == "usemtl") {
                file >> current_material;
            }
        }
    }

    void AddTriangle(Triangle triangle) { triangles_.push_back(triangle); }

    void TraceRay(Ray &ray) {
        float distance = INFINITY;
        int triangles_count = triangles_.size();
        int i = 0, current_triangle = -1;
        for (i = 0; i < triangles_count; ++i) {
            if (triangles_[i].Intersect(ray, distance)) {
                current_triangle = i;
            }
        }
        if (current_triangle > -1) {
            vec4 drop_point = ray.GetBegin() + ray.GetDirection() * distance;
            vec4 normal = triangles_[current_triangle].GetNormal();
            triangles_[current_triangle].GetMaterial().Process(ray, drop_point,
                                                               normal);
        } else {
            if (Config::get().skybox != "") {
                float theta = acos(ray.GetDirection().y) / pi,
                        phi =
                        atan2(ray.GetDirection().z, -ray.GetDirection().x) / pi / 2 +
                        0.5f;

                float x = (phi * skybox_.width()), y = (theta * skybox_.height());
                int x1 = (int)x, y1 = (int)y;
                int x2 = (x1 + 1) % skybox_.width(), y2 = (y1 + 1) % skybox_.height();

                auto _color1 = skybox_.get_pixel(x1, y1);
                auto _color2 = skybox_.get_pixel(x2, y1);
                auto _color3 = skybox_.get_pixel(x1, y2);
                auto _color4 = skybox_.get_pixel(x2, y2);

                vec3 color1 =
                        vec3((float)_color1.red, (float)_color1.green, (float)_color1.blue);
                vec3 color2 =
                        vec3((float)_color2.red, (float)_color2.green, (float)_color2.blue);
                vec3 color3 =
                        vec3((float)_color3.red, (float)_color3.green, (float)_color3.blue);
                vec3 color4 =
                        vec3((float)_color4.red, (float)_color4.green, (float)_color4.blue);

                vec3 color12 = mix(color1, color2, 1 - x + x1);
                vec3 color34 = mix(color3, color4, 1 - x + x1);

                vec3 color = mix(color12, color34, 1 - y + y1) / 256.f;

                color_map_[ray.GetCoords().x][ray.GetCoords().y] += color;
                color2_map_[ray.GetCoords().x][ray.GetCoords().y] += color * color;
                ++samples_count_[ray.GetCoords().x][ray.GetCoords().y];
            }
            ray.MakeInvalid();
        }
    }
};

int main(int argc, char *argv[]) {
    chrono::milliseconds start_time = chrono::duration_cast<chrono::milliseconds>(
            chrono::system_clock::now().time_since_epoch());

    Config::get().set_config(argc, argv);
    static default_random_engine generator(time(0));
    static uniform_real_distribution<> distribution(-0.5f, 0.5f);

    vector<vector<vec3> > color_map(Config::get().width,
                                    vector<vec3>(Config::get().height, vec3(0)));
    vector<vector<vec3> > color2_map(Config::get().width,
                                     vector<vec3>(Config::get().height, vec3(0)));
    vector<vector<int> > samples_count(Config::get().width,
                                       vector<int>(Config::get().height, 0));

    Scene scene = Scene(color_map, color2_map, samples_count);
    scene.LoadModel(Config::get().model_path + Config::get().model_name);

    bitmap_image image(Config::get().width, Config::get().height);

    image.clear();
    for (int i = 0; i < Config::get().rays_per_pixel; ++i) {
        for (int y = 0; y < Config::get().height; ++y) {
            #pragma omp parallel for num_threads(4)
            for (int x = 0; x < Config::get().width; ++x) {
                float sample_count = static_cast<float>(samples_count[x][y]);
                if (i > 10 && sample_count) {
                    vec3 D = (color2_map[x][y] / sample_count -
                              (color_map[x][y] / sample_count) *
                              (color_map[x][y] / sample_count));
                    if ((i % 4) && D.r < Config::get().error &&
                        D.g < Config::get().error && D.b < Config::get().error) {
                        continue;
                    }
                }
                vec4 direction = vec4(
                        (x + distribution(generator)) / Config::get().width - 0.5f,
                        -(y + distribution(generator)) / Config::get().height + 0.5f, 1, 0);
                Ray ray = Ray(vec4(0, 0, -20, 1), direction, 0, ivec2(x, y));
                while (ray.IsValid()) {
                    scene.TraceRay(ray);
                }
                if (i % Config::get().update == 0) {
                    if (samples_count[x][y]) {
                        vec3 c =
                                pow(color_map[x][y] / static_cast<float>(samples_count[x][y]),
                                    vec3(Config::get().gamma_correction)) *
                                255.0f;
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
            if (!samples_count[x][y]) {
                continue;
            }
            color_map[x][y] =
                    pow(color_map[x][y] / static_cast<float>(samples_count[x][y]),
                        vec3(1 / 2.2)) *
                    255.0f;
        }
    }

    if (Config::get().gauss) {
        color_map = GaussBlur(color_map, Config::get().gauss);
    }
    if (Config::get().median) {
        color_map = MedianFilter(color_map, Config::get().median);
    }
    for (int y = 0; y < Config::get().height; ++y) {
        for (int x = 0; x < Config::get().width; ++x) {
            if (!samples_count[x][y]) {
                continue;
            }
            vec3 color = color_map[x][y];
            image.set_pixel(x, y, color.r, color.g, color.b);
        }
    }

    chrono::milliseconds end_time = chrono::duration_cast<chrono::milliseconds>(
            chrono::system_clock::now().time_since_epoch());
    time_t t = time(0);
    struct tm *now = localtime(&t);
    string date = to_string(now->tm_year + 1900) + '-' +
                  to_string(now->tm_mon + 1) + '-' + to_string(now->tm_mday) +
                  '-' + to_string(now->tm_hour) + '-' + to_string(now->tm_min) +
                  '-' + to_string(now->tm_sec) + "  " +
                  to_string((end_time - start_time).count()) + "   " +
                  to_string(Config::get().rays_per_pixel);
    image.save_image(date + ".bmp");
    image.save_image("result.bmp");
}