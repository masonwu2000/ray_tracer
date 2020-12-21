// header files
#include "src/rtweekend.h"
#include "src/color.h"
#include "src/hittable_list.h"
#include "src/sphere.h"
#include "src/camera.h"
#include "src/material.h"

#include <iostream>

hittable_list random_scene() {
    hittable_list world;

    shared_ptr<lambertian> ground = make_shared<lambertian>(color(0.5, 0.5, 0.5));
    world.add(make_shared<sphere>(point3(0.0, -100.0, 0.0), 100.0, ground));

    for (int a = -5; a < 5; a++) {
        for (int b = -5; b < 5; b++) {
            double choose_mat = random_double();
            point3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());

            if ((center - point3(4.0, 0.2, 0.0)).length() > 0.9) {
                shared_ptr<material> sphere_material;

                if (choose_mat < 0.8) {
                    vec3 albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else if (choose_mat < 0.95) {
                    vec3 albedo = color::random(0.5, 1.0);
                    double fuzz = random_double(0.0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else {
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    shared_ptr<dielectric> material1 = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0.0, 1.0, 0.0), 1.0, material1));

    shared_ptr<lambertian> material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4.0, 1.0, 0.0), 1.0, material2));

    shared_ptr<metal> material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4.0, 1.0, 0.0), 1.0, material3));

    return world;
}

double hit_sphere(const point3& center, double radius, const ray& r) {
    vec3 oc = r.origin() - center;
    double a = dot(r.direction(), r.direction());
    double half_b = dot(r.direction(), oc);
    double c = dot(oc, oc) - radius * radius;
    double discriminant = half_b * half_b - a * c;
    if (discriminant < 0) {
        return -1.0;
    }
    return (-half_b - sqrt(discriminant)) / a;
}

color ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;
    if (depth <= 0) {
        return color(0.0, 0.0, 0.0);
    }
    if (world.hit(r, 0.001, infinity, rec)) {
        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
            return attenuation * ray_color(scattered, world, depth - 1);
        }
        return color(0.0, 0.0, 0.0);
    }
    vec3 unit_direction = unit_vector(r.direction());
    double t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
}

int main() {

    // Image size
    const double aspect_ratio = 3.0 / 2.0;
    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 10;
    const int max_depth = 50;

    // World
    hittable_list world = random_scene();

    point3 lookfrom(13.0, 2.0, 3.0);
    point3 lookat(0.0, 0.0, 0.0);
    vec3 vup(0.0, 1.0, 0.0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.1;

    camera cam(lookfrom, lookat, vup, 20.0, aspect_ratio, aperture, dist_to_focus);

    // Render
    
    cout << "P3\n" << image_width << ' ' << image_height << "\n255" << endl;

    for (int i = image_height - 1; i >= 0 ; i--) {
        
        cerr << "Lines remaining: " << i << endl << flush;

        for (int j = 0; j < image_width; j++) {
            color pixel_color(0.0, 0.0, 0.0);
            for (int s = 0; s < samples_per_pixel; s++) {
                double u = (j + random_double()) / (image_width - 1);
                double v = (i + random_double()) / (image_height - 1);

                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world, max_depth);
            }
            write_color(cout, pixel_color, samples_per_pixel);
        }
    }

    cerr << "Complete" << endl;
}
