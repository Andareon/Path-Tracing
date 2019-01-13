//
// Created by alex on 13.01.2019.
//

#include "scene.h"

using namespace std;
using namespace glm;

static vector<string> Split(const string &str, char sep) {
    vector<string> ans;
    int begin, end;
    for (begin = 0, end = str.find(sep); end != string::npos; begin = end + 1, end = str.find(sep, begin)) {
        ans.push_back(str.substr(begin, end - begin));
    }
    ans.push_back(str.substr(begin));
    return ans;
}

Scene::Scene(vector<vector<vec3>> &color_map, vector<vector<vec3>> &color2_map, vector<vector<int> > &samples_count)
: color_map_(color_map),
color2_map_(color2_map),
samples_count_(samples_count) {
    if (!Config::get().skybox.empty()) {
        skybox_ = bitmap_image(Config::get().skybox);
    }
};


void Scene::LoadModel(string path) {
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
            MaterialCharacteristics characteristics = {0};
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

void Scene::AddTriangle(Triangle triangle) { triangles_.push_back(triangle); }

void Scene::TraceRay(Ray &ray) {
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
        triangles_[current_triangle].GetMaterial().Process(ray, drop_point, normal);
    } else {
        if (!Config::get().skybox.empty()) {
            const float theta = acos(ray.GetDirection().y) / pi;
            const float phi = atan2(ray.GetDirection().z, -ray.GetDirection().x) / pi / 2 + 0.5f;

            const float x = (phi * skybox_.width()), y = (theta * skybox_.height());
            const auto x1 = static_cast<const unsigned int>(x);
            const auto y1 = static_cast<const unsigned int>(y);
            const unsigned x2 = (x1 + 1) % skybox_.width();
            const unsigned y2 = (y1 + 1) % skybox_.height();

            const auto _color1 = skybox_.get_pixel(x1, y1);
            const auto _color2 = skybox_.get_pixel(x2, y1);
            const auto _color3 = skybox_.get_pixel(x1, y2);
            const auto _color4 = skybox_.get_pixel(x2, y2);

            const vec3 color1 = vec3(_color1.red, _color1.green, _color1.blue);
            const vec3 color2 = vec3(_color2.red, _color2.green, _color2.blue);
            const vec3 color3 = vec3(_color3.red, _color3.green, _color3.blue);
            const vec3 color4 = vec3(_color4.red, _color4.green, _color4.blue);

            const vec3 color12 = mix(color1, color2, 1 - x + x1);
            const vec3 color34 = mix(color3, color4, 1 - x + x1);

            const vec3 color = mix(color12, color34, 1 - y + y1) / 256.f;

            color_map_[ray.GetCoords().x][ray.GetCoords().y] += color;
            color2_map_[ray.GetCoords().x][ray.GetCoords().y] += color * color;
            ++samples_count_[ray.GetCoords().x][ray.GetCoords().y];
        }
        ray.MakeInvalid();
    }
}