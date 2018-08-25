#include <fstream>
#include <cmath>
#include <cstdlib>
#include "bitmap_image.hpp"

using namespace std;

struct Vec {
    float x, y, z;
    Vec(){x=y=z=0;}
    Vec(float a, float b, float c) {
        x=a, y=b, z=c;
    }
    Vec operator- (Vec v) const {
        return Vec(x - v.x, y - v.y, z - v.z);
    }

    Vec operator* (float b) const {
        return Vec(x*b, y*b, z*b);
    }

    Vec operator/ (float b) const {
        return Vec(x/b, y/b, z/b);
    }
    Vec operator+ (Vec v) const {
        return Vec(x + v.x, y + v.y, z + v.z);
    }

    Vec normalize() const {
        float mg = sqrt(x*x+y*y+z*z);
        return Vec(x/mg, y/mg, z/mg);
    }

};

float dot(Vec v, Vec b) {
    return v.x*b.x+v.y*b.y+v.z*b.z;
}

struct Ray {
    Vec o;
    Vec d;
    Ray(Vec i, Vec j) {o=i, d=j;}
};

struct Color {
    float r, g, b;
    Color(){r=g=b=0;}
    Color(float i, float j, float k){r=i,g=j,b=k;}
    Color operator* (float d) const {
        return Color(r*d, g*d, b*d);
    }
    Color operator+ (Color c) const {
        return Color(r+c.r, g+c.g, b+c.b);
    }
};

struct Sphere {
    Vec c;
    float r;
    Color col;
    Sphere(Vec i, float j, Color k){c=i, r=j, col=k;}

    Vec getNormal(Vec pi) const {
        return (pi - c) / r;
    }

    bool intersect(const Ray &ray, float &t) const {
        Vec o = ray.o;
        Vec d = ray.d.normalize();
        Vec oc = o-c;
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

float saturate(float x) {
    return std::max(std::min(255.0f, x), 0.0f);
}

Color saturate(Color x) {
    return Color(saturate(x.r), saturate(x.g), saturate(x.b));
}

Color traceRay(const Ray& ray) {
    static const Color white (255, 255, 255);
    static const Color red (255, 0, 0);
    static const Color blue (128, 0, 255);
    static const Color violet (255, 0, 255);
    static const Color yellow (255, 255, 0);
    static const Color black (0, 0, 0);
    
    static const Sphere spheres[3] = {Sphere(Vec(-11,0,20),5, red), Sphere(Vec(0,0,35),5, blue), Sphere(Vec(11,0,25),5, violet)};
    static const Sphere light (Vec(0,20,15),1, white);
    
    Color result = black;
    
    float t = 20000000;
    for (const auto &sphere: spheres) {
        if (sphere.intersect(ray, t)) {
            Vec pi = ray.o + ray.d * t;

            Vec L = light.c - pi;
            Vec N = sphere.getNormal(pi);
            float dt = std::max(dot(L.normalize(), N.normalize()), 0.0f);



            Vec NorL = L.normalize();
            Vec NorN = N.normalize();
            Vec ans = NorL * -1 + NorN * 2 * dot(NorN, NorL);
            float cs = std::max(0.0f, (dot(ans, pi * -1) / t));







            dt += std::min(0.05f, (1 - dt));

            result = sphere.col * dt + light.col * std::pow(cs, 32);


        }
    }
    return saturate(result);
}

float randFloat(const float min, const float max) {
    return (static_cast<float>(rand()) / RAND_MAX) * (max - min) + min;
}


int main() {
    int W = 500;
    int H = 500;
    int ray_col = 200;

    bitmap_image image(H, W);
    image.clear();

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            float r = 0, g = 0, b = 0;
            for (int i = 0; i < ray_col; ++i) {
                Ray ray(Vec(0, 0, 0), Vec(static_cast<float>(x - W / 2 + randFloat(-0.5f, 0.5f)) / W * 2,
                                          -static_cast<float>(y - H / 2 + randFloat(-0.5f, 0.5f)) / H * 2,
                                          1).normalize());
                Color c = traceRay(ray);
                r += c.r;
                g += c.g;
                b += c.b;
            }
                image.set_pixel(x, y, r / ray_col, g / ray_col, b / ray_col);

        }
    }

//    W /= 4;
//    H /= 4;
//    bitmap_image image2(H, W);
//    image2.clear();
//    for (int y = 0; y < H; ++y) {
//        for (int x = 0; x < W; ++x)  {
//            int r = 0, g = 0, b = 0;
//            for (int dx = 0; dx < 4; ++dx) {
//                for (int dy = 0; dy < 4; ++dy) {
//                    auto c = image.get_pixel(x * 4 + dx, y * 4 + dy);
//                    r += c.red;
//                    g += c.green;
//                    b += c.blue;
//                }
//            }
//            image2.set_pixel(x, y, r / 16, g / 16, b / 16);
//        }
//    }


    image.save_image("1.bmp");
}