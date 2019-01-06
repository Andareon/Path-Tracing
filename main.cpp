#include <chrono>

#include "bitmap_image.hpp"

#include "material.h"

using namespace std;
using namespace glm;

bool planeintersect(Ray &ray, float &t, glm::vec4 plane) {
    glm::vec4 o = ray.getBegin();
    glm::vec4 d = ray.getDir();
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

float square(glm::vec4 A, glm::vec4 B, glm::vec4 C) {
    glm::vec3 a = B - A;
    glm::vec3 b = C - A;
    glm::vec3 c = cross(a, b);
    return length(c) / 2;
}

glm::vec4 cross(glm::vec4 a, glm::vec4 b) {
    return glm::vec4(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x, 0);
}

vector<string> split(string str, char sep) {
    vector<string> ans;
    int begin, end;
    for (begin = 0, end = str.find(sep); end != std::string::npos; begin = end + 1, end = str.find(sep, begin)) {
        ans.push_back(str.substr(begin, end - begin));
    }
    ans.push_back(str.substr(begin));
    return ans;
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

class Triangle {
private:
    glm::vec4 plane;
    Material material;
    std::array<glm::vec4, 3> vertices;

public:
    Triangle(const std::array<glm::vec4, 3> &a, Material m): vertices(a), material(m) {
        glm::vec4 e1 = vertices[1] - vertices[0];
        glm::vec4 e2 = vertices[2] - vertices[0];
        plane = normalize(cross(e1, e2));
        plane.w = -dot(plane, vertices[0]);
    }

    glm::vec4 getNormal() const {
        return glm::vec4(plane.x, plane.y, plane.z, 0);
    }

    void setNormal(glm::vec3 N) {
        N = normalize(N);
        plane.x = N.x;
        plane.y = N.y;
        plane.z = N.z;
        plane.w = 0;
        plane.w = -dot(plane, vertices[0]);
    }

    Material getMaterial() const {
        return material;
    }

    bool intersect(Ray &ray, float &t) const {
        float newT = INFINITY;
        if (!planeintersect(ray, newT, plane)) {
            return false;
        }

        if (newT >= t) {
            return false;
        }

        glm::vec4 drop_point = ray.getBegin() + ray.getDir() * newT;
        float full_square = square(vertices[0], vertices[1], vertices[2]);
        float small_square_1 = square(drop_point, vertices[1], vertices[2]);
        float small_square_2 = square(vertices[0], drop_point, vertices[2]);
        float small_square_3 = square(vertices[0], vertices[1], drop_point);
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
    vector<vector<vec3> > &ColorMap;
    vector<vector<vec3> > &Color2Map;
    vector<vector<int> > &SamplesCount;
    vector<Material> materials;
    bitmap_image skybox;
public:
    Scene(vector<vector<vec3> > &CM, vector<vector<vec3> > &C2M, vector<vector<int> > &SC) :ColorMap(CM), Color2Map(C2M), SamplesCount(SC){
        if (Config::get().skybox != "") {
            skybox = bitmap_image(Config::get().skybox);
        }
    };
    void LoadModel(string path) {
        vector<vec4> temp_vertices;
        vector<vec2> temp_texture_coords;
        vector<vec3> temp_normals;
        ifstream file(path);
        if (file.eof()) {
            printf("Impossible to open the file !\n");
            return;
        }
        int curmat = 0;
        while (1) {
            if (file.eof())
                break;
            string first;
            file >> first;
            if (first == "mtllib") {
                string mtlname;
                file >> mtlname;
                ifstream mtlfile(Config::get().model_path + mtlname);

                string mtlfirst = "1";
                Material_characteristics characteristics;
                while(1) {
                    if (mtlfile.eof()) {
                        break;
                    }
                    while (!mtlfile.eof() && mtlfirst != "newmtl") {
                        mtlfile >> mtlfirst;
                    }
                    mtlfile >> mtlfirst;
                    while (!mtlfile.eof() && mtlfirst != "newmtl") {
                        if (mtlfirst == "Kd") {
                            mtlfile >> characteristics.Kd.x >> characteristics.Kd.y >> characteristics.Kd.z;
                        } else if (mtlfirst == "Ke") {
                            mtlfile >> characteristics.Ke.x >> characteristics.Ke.y >> characteristics.Ke.z;
                        }
                        mtlfile >> mtlfirst;
                    }
                    materials.emplace_back(Factory(characteristics, ColorMap, Color2Map, SamplesCount));
                }

            } else if (first == "v") {
                glm::vec4 vertex = vec4(1);
                file >> vertex.x >> vertex.y >> vertex.z;
                temp_vertices.push_back(vertex);
            } else if (first == "vt") {
                glm::vec2 uv;
                file >> uv.x >> uv.y;
                temp_texture_coords.push_back(uv);
            } else if (first == "vn") {
                glm::vec3 normal;
                file >> normal.x >> normal.y >> normal.z;
                temp_normals.push_back(normal);
            } else if (first == "f") {
                int vertex_index[3];
                int texture_coords_index[3];
                int normal_index[3];
                for (int i = 0; i < 3; ++i) {
                    string cr;
                    file >> cr;
                    vector<string> cur = split(cr, '/');
                    cur.resize(3);
                    vertex_index[i] = atoi(cur[0].c_str()) - 1;
                    texture_coords_index[i] = atoi(cur[1].c_str()) - 1;
                    normal_index[i] = atoi(cur[2].c_str()) - 1;
                }
                triangles.emplace_back(array<vec4, 3>{temp_vertices[vertex_index[0]], temp_vertices[vertex_index[1]], temp_vertices[vertex_index[2]]}, materials[curmat]);
                if (normal_index[0] >= 0) {
                    triangles.back().setNormal(temp_normals[normal_index[0]]);
                }
            } else if (first == "usemtl") {
                file >> curmat;
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
            vec4 drop_point = ray.getBegin() + ray.getDir() * t;
            vec4 N = triangles[cur].getNormal();
            triangles[cur].getMaterial().process(ray, drop_point, N);
        } else {
            if (Config::get().skybox != "") {
                float theta = acos(ray.getDir().y) / PI,
                        phi = atan2(ray.getDir().z, -ray.getDir().x) / PI / 2 + 0.5f;

                float x = (phi * skybox.width()), y = (theta * skybox.height());
                int x1 = (int) x, y1 = (int) y;
                int x2 = (x1 + 1) % skybox.width(), y2 = (y1 + 1) % skybox.height();

                auto col1 = skybox.get_pixel(x1, y1);
                auto col2 = skybox.get_pixel(x2, y1);
                auto col3 = skybox.get_pixel(x1, y2);
                auto col4 = skybox.get_pixel(x2, y2);

                vec3 cl1 = vec3((float) col1.red, (float) col1.green, (float) col1.blue);
                vec3 cl2 = vec3((float) col2.red, (float) col2.green, (float) col2.blue);
                vec3 cl3 = vec3((float) col3.red, (float) col3.green, (float) col3.blue);
                vec3 cl4 = vec3((float) col4.red, (float) col4.green, (float) col4.blue);

                vec3 cl12 = mix(cl1, cl2, 1 - x + x1);
                vec3 cl34 = mix(cl3, cl4, 1 - x + x1);

                vec3 cl = mix(cl12, cl34, 1 - y + y1) / 256.f;

                ColorMap[ray.getCoords().x][ray.getCoords().y] += cl;
                Color2Map[ray.getCoords().x][ray.getCoords().y] += cl * cl;
                ++SamplesCount[ray.getCoords().x][ray.getCoords().y];
            }
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

    Scene scene = Scene(ColorMap, Color2Map, SamplesCount);
    scene.LoadModel(Config::get().model_path + Config::get().model_name);


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
                vec4 dir = vec4((x + distribution(generator)) / Config::get().width - 0.5f,
                                -(y + distribution(generator)) / Config::get().height + 0.5f,
                                1, 0);
                Ray ray = Ray(vec4(0, 0, -20, 1), dir, 0, ivec2(x, y));
                while (ray.is_valid()) {
                    scene.traceRay(ray);
                }
                if (i % Config::get().update == 0) {
                    if (SamplesCount[x][y]) {
                        vec3 c = pow(ColorMap[x][y] / static_cast<float>(SamplesCount[x][y]), vec3(Config::get().gamma_kor)) * 255.0f;
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


    if (Config::get().gauss) {
        ColorMap = gauss_blur(ColorMap, Config::get().gauss);
    }
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