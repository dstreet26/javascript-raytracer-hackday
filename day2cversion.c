

// todo: world with random spheres and materials
// todo: remove intermediate steps, use more types, see if allocating memory is needed

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

enum MaterialType
{
    LAMBERTIAN,
    METAL,
    GLASS
};

// todo: would it help to make things constant?
struct Material
{
    enum MaterialType type;
    float albedo[3];
    float ir;
    float fuzz;
};
struct Sphere
{
    float center[3];
    float radius;
    struct Material material;
};
struct Ray
{
    float origin[3];
    float direction[3];
};
struct Hit
{
    float t;
    float p[3];
    int frontFace;
    float normal[3];
    struct Material material;
};

struct Scatter
{
    struct Ray scatteredRay;
    float attenuation[3];
};
struct Rec
{
    struct Hit hit;
    struct Material material;
};
void subtractVectors(const float a[3], const float b[3], float result[3])
{
    result[0] = a[0] - b[0];
    result[1] = a[1] - b[1];
    result[2] = a[2] - b[2];
}
void addVectors(float a[3], float b[3], float result[3])
{
    result[0] = a[0] + b[0];
    result[1] = a[1] + b[1];
    result[2] = a[2] + b[2];
}
void multiplyVector(float a[3], float b, float result[3])
{
    result[0] = a[0] * b;
    result[1] = a[1] * b;
    result[2] = a[2] * b;
}
void multiplyVectors(float a[3], float b[3], float result[3])
{
    result[0] = a[0] * b[0];
    result[1] = a[1] * b[1];
    result[2] = a[2] * b[2];
}

float vectorLength(float vector[3])
{
    return sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
}
void unitVector(float vector[3], float result[3])
{
    float length = vectorLength(vector);
    multiplyVector(vector, 1 / length, result);
}

float dotVectors(const float a[3], const float b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
void crossVectors(const float a[3], const float b[3], float result[3])
{
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = -(a[0] * b[2] - a[2] * b[0]);
    result[2] = a[0] * b[1] - a[1] * b[0];
}

float randomFloat(float min, float max)
{
    return (float)rand() / RAND_MAX * (max - min) + min;
}
void randomVector3(float min, float max, float result[3])
{
    result[0] = randomFloat(min, max);
    result[1] = randomFloat(min, max);
    result[2] = randomFloat(min, max);
    return;
}
float lengthSquared(float vector[3])
{
    return vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2];
}
void random_in_unit_sphere(float result[3])
{
    while (1)
    {
        randomVector3(-1, 1, result);
        if (lengthSquared(result) < 1)
        {
            return;
        }
    }
}
void random_unit_vector(float result[3])
{
    float p[3];
    random_in_unit_sphere(p);
    unitVector(p, result);
    return;
}
void reflect(float v[3], float n[3], float result[3])
{
    float temp[3];
    multiplyVector(n, 2 * dotVectors(v, n), temp);
    subtractVectors(v, temp, result);
    return;
}


void refract(float uv[3], float n[3], float etai_over_etat, float result[3])
{
    float intermediatestep1[3];
    multiplyVector(uv, -1, intermediatestep1);
    float cos_theta = fmin(dotVectors(intermediatestep1, n), 1.0);
    float r_out_perp[3];
    float intermediatestep2[3];
    multiplyVector(n, cos_theta, intermediatestep2);
    addVectors(uv, intermediatestep2, r_out_perp);
    multiplyVector(r_out_perp, etai_over_etat, r_out_perp);
    float r_out_parallel[3];
    multiplyVector(n, -sqrt(fabs(1.0 - lengthSquared(r_out_perp))), r_out_parallel);
    addVectors(r_out_perp, r_out_parallel, result);
    return;
}
float reflectance(float cosine, float ref_idx)
{
    float r0 = (1 - ref_idx) / (1 + ref_idx);
    r0 = r0 * r0;
    return r0 + (1 - r0) * pow(1 - cosine, 5);
}

struct Ray getRay(int i, int j, float pixel00Location[3], float pixelDeltaU[3], float pixelDeltaV[3], const float cameraCenter[3])
{
    float pixelCenter[3];
    float intermediatestep1[3];
    float intermediatestep2[3];
    float intermediatestep3[3];
    multiplyVector(pixelDeltaU, i, intermediatestep1);
    multiplyVector(pixelDeltaV, j, intermediatestep2);
    addVectors(intermediatestep1, intermediatestep2, intermediatestep3);
    addVectors(pixel00Location, intermediatestep3, pixelCenter);
    float direction[3];
    subtractVectors(pixelCenter, cameraCenter, direction);
    struct Ray ray = {
        .origin = {cameraCenter[0], cameraCenter[1], cameraCenter[2]},
        .direction = {direction[0], direction[1], direction[2]},
    };
    return ray;
}

int surrounds(float x, float min, float max)
{
    return min <= x && x <= max;
}
struct Hit hitSphere(float center[3], float radius, struct Ray ray, float min, float max)

{
    float oc[3];
    subtractVectors(ray.origin, center, oc);
    float a = lengthSquared(ray.direction);

    float halfB = dotVectors(oc, ray.direction);
    float c = lengthSquared(oc) - radius * radius;
    float discriminant = halfB * halfB - a * c;
    if (discriminant < 0)
    {
        struct Hit hit = {
            .t = -1,
        };
        return hit;
    }
    else
    {
        float sqrtd = sqrt(discriminant);
        float root = (-halfB - sqrtd) / a;
        if (!surrounds(root, min, max))
        {
            root = (-halfB + sqrtd) / a;
            if (!surrounds(root, min, max))
            {
                struct Hit hit = {
                    .t = -1,
                };
                return hit;
            }
        }
        float p[3];
        multiplyVector(ray.direction, root, p);
        addVectors(ray.origin, p, p);
        float outwardNormal[3];
        subtractVectors(p, center, outwardNormal);
        multiplyVector(outwardNormal, 1 / radius, outwardNormal);
        int frontFace = dotVectors(ray.direction, outwardNormal) < 0;
        float normal[3];
        if (frontFace)
        {
            normal[0] = outwardNormal[0];
            normal[1] = outwardNormal[1];
            normal[2] = outwardNormal[2];
        }
        else
        {
            normal[0] = -outwardNormal[0];
            normal[1] = -outwardNormal[1];
            normal[2] = -outwardNormal[2];
        }
        struct Hit hit = {
            .t = root,
            .p = {p[0], p[1], p[2]},
            .frontFace = frontFace,
            .normal = {normal[0], normal[1], normal[2]},
        };
        return hit;
    }
}

struct Scatter scatter(struct Ray ray, struct Rec rec)
{
    struct Scatter scatter;
    if (rec.material.type == LAMBERTIAN)
    {
        scatter.scatteredRay.origin[0] = rec.hit.p[0];
        scatter.scatteredRay.origin[1] = rec.hit.p[1];
        scatter.scatteredRay.origin[2] = rec.hit.p[2];
        float randomUnitVector[3];
        random_unit_vector(randomUnitVector);
        scatter.scatteredRay.direction[0] = rec.hit.normal[0] + randomUnitVector[0];
        scatter.scatteredRay.direction[1] = rec.hit.normal[1] + randomUnitVector[1];
        scatter.scatteredRay.direction[2] = rec.hit.normal[2] + randomUnitVector[2];
        scatter.attenuation[0] = rec.material.albedo[0];
        scatter.attenuation[1] = rec.material.albedo[1];
        scatter.attenuation[2] = rec.material.albedo[2];
    }
    else if (rec.material.type == METAL)
    {
        float reflected[3];
        reflect(ray.direction, rec.hit.normal, reflected);
        float fuzzed[3];
        float randomUnitVector[3];
        random_unit_vector(randomUnitVector);
        multiplyVector(randomUnitVector, rec.material.fuzz, fuzzed);
        float scattered[3];
        addVectors(reflected, fuzzed, scattered);
        if (dotVectors(scattered, rec.hit.normal) > 0)
        {
            scatter.scatteredRay.origin[0] = rec.hit.p[0];
            scatter.scatteredRay.origin[1] = rec.hit.p[1];
            scatter.scatteredRay.origin[2] = rec.hit.p[2];
            scatter.scatteredRay.direction[0] = scattered[0];
            scatter.scatteredRay.direction[1] = scattered[1];
            scatter.scatteredRay.direction[2] = scattered[2];
            scatter.attenuation[0] = rec.material.albedo[0];
            scatter.attenuation[1] = rec.material.albedo[1];
            scatter.attenuation[2] = rec.material.albedo[2];
        }
    }
    else if (rec.material.type == GLASS)
    {
        float refractionRatio = rec.hit.frontFace ? 1 / rec.material.ir : rec.material.ir;
        float unitDirection[3];
        unitVector(ray.direction, unitDirection);
        float intermediatestep1[3];
        multiplyVector(unitDirection, -1, intermediatestep1);
        float cosTheta = fmin(dotVectors(intermediatestep1, rec.hit.normal), 1);

        float sinTheta = sqrt(1.0 - cosTheta * cosTheta);
        float cannotRefract = refractionRatio * sinTheta > 1.0;
        float direction[3];
        if (cannotRefract || reflectance(cosTheta, refractionRatio) > 0.5)
        {
            reflect(unitDirection, rec.hit.normal, direction);
        }
        else
        {
            refract(unitDirection, rec.hit.normal, refractionRatio, direction);
        }
        scatter.scatteredRay.origin[0] = rec.hit.p[0];
        scatter.scatteredRay.origin[1] = rec.hit.p[1];
        scatter.scatteredRay.origin[2] = rec.hit.p[2];
        scatter.scatteredRay.direction[0] = direction[0];
        scatter.scatteredRay.direction[1] = direction[1];
        scatter.scatteredRay.direction[2] = direction[2];
        scatter.attenuation[0] = rec.material.albedo[0];
        scatter.attenuation[1] = rec.material.albedo[1];
        scatter.attenuation[2] = rec.material.albedo[2];
    }
    return scatter;
}

void rayColor(struct Ray ray, float depth, float result[3], struct Sphere world[], int worldLength)

{

    if (depth <= 0)
    {
        result[0] = 0;
        result[1] = 0;
        result[2] = 0;
        return;
    }
    float min1 = 0.001;
    float max1 = 100000;

    struct Rec rec;
    rec.hit.t = -1;
    for (int k = 0; k < worldLength; k++)
    {
        struct Sphere element = world[k];
        struct Hit hit = hitSphere(element.center, element.radius, ray, min1, max1);
        if (hit.t != -1)
        {
            max1 = hit.t;
            rec.hit = hit;

            rec.material = element.material;
        }
    }
    if (rec.hit.t != -1)
    {

        struct Scatter scattered = scatter(ray, rec);
        float tempcolor[3];
        rayColor(scattered.scatteredRay, depth - 1, tempcolor, world, worldLength);
        multiplyVectors(tempcolor, scattered.attenuation, result);
        return;
    }
    float unitDirection[3];
    unitVector(ray.direction, unitDirection);
    float a = 0.5 * (unitDirection[1] + 1.0);
    result[0] = (1.0 - a) * 1.0 + a * 0.5;
    result[1] = (1.0 - a) * 1.0 + a * 0.7;
    result[2] = (1.0 - a) * 1.0 + a * 1.0;
    return;
}
void writeColor(FILE *f, float color[3], int samplesPerPixel)
{
    float scale = 1.0 / samplesPerPixel;
    float scaledColor[3];
    multiplyVector(color, scale, scaledColor);
    scaledColor[0] = sqrt(scaledColor[0]);
    scaledColor[1] = sqrt(scaledColor[1]);
    scaledColor[2] = sqrt(scaledColor[2]);

    float colorScaled[3];
    multiplyVector(scaledColor, 255, colorScaled);
    colorScaled[0] = round(colorScaled[0]);
    colorScaled[1] = round(colorScaled[1]);
    colorScaled[2] = round(colorScaled[2]);
    fprintf(f, "%f %f %f\n", colorScaled[0], colorScaled[1], colorScaled[2]);
}
int main()
{
    const int imageWidth = 1000;
    const int fov = 20;
    const int focusDist = 20;
    const int samplesPerPixel = 10;
    const int maxDepth = 5;
    const int maxVal = 255;
    const float aspectRatio = 16.0 / 9.0;
    const int imageHeight = floor(imageWidth / aspectRatio);
    const float viewportHeight = 2 * tan((fov * (M_PI / 180)) / 2) * focusDist;
    const float viewportWidth = viewportHeight * ((float)imageWidth / imageHeight);
    const float cameraCenter[3] = {13, 2, 3};
    const float lookAt[3] = {0, 0, 0};
    const float vUp[3] = {0, 1, 0};
    float w[3];
    subtractVectors(cameraCenter, lookAt, w);
    unitVector(w, w);
    float u[3];
    crossVectors(vUp, w, u);
    unitVector(u, u);
    float v[3];
    crossVectors(w, u, v);
    float viewportU[3];
    multiplyVector(u, viewportWidth, viewportU);
    float viewportV[3];
    multiplyVector(v, -viewportHeight, viewportV);
    float pixelDeltaU[3];
    multiplyVector(viewportU, 1.0 / imageWidth, pixelDeltaU);
    float pixelDeltaV[3];
    multiplyVector(viewportV, 1.0 / imageHeight, pixelDeltaV);
    float viewportUpperLeft[3];
    float intermediatestep1[3];
    float intermediatestep2[3];
    multiplyVector(w, focusDist, intermediatestep1);
    subtractVectors(cameraCenter, intermediatestep1, intermediatestep2);
    multiplyVector(viewportU, 0.5, intermediatestep1);
    subtractVectors(intermediatestep2, intermediatestep1, intermediatestep2);
    multiplyVector(viewportV, 0.5, intermediatestep1);
    subtractVectors(intermediatestep2, intermediatestep1, viewportUpperLeft);

    float pixel00Location[3];
    float intermediatestep3[3];
    float intermediatestep4[3];
    addVectors(pixelDeltaU, pixelDeltaV, intermediatestep3);
    multiplyVector(intermediatestep3, 0.5, intermediatestep4);
    addVectors(viewportUpperLeft, intermediatestep4, pixel00Location);

    // todo: would it help to make things constant?
    struct Sphere groundSphere = {
        .center = {0, -1000, 0},
        .radius = 1000,
        .material = {
            .type = LAMBERTIAN,
            .albedo = {0.5, 0.5, 0.5},
        },
    };
    struct Sphere referenceSphere1 = {
        .center = {0, 1, 0},
        .radius = 1,
        .material = {
            .type = GLASS,
            .albedo = {1, 1, 1},
            .ir = 1.5,
        },
    };
    struct Sphere referenceSphere2 = {
        .center = {-4, 1, 0},
        .radius = 1,
        .material = {
            .type = LAMBERTIAN,
            .albedo = {0.4, 0.2, 0.1},
        },
    };

    struct Sphere referenceSphere3 = {
        .center = {4, 1, 0},
        .radius = 1,
        .material = {
            .type = METAL,
            .albedo = {0.7, 0.6, 0.5},
            .fuzz = 0,
        },
    };

    const int minA = -11;
    const int maxA = 11;
    const int minB = -11;
    const int maxB = 11;
    const int count = (maxA - minA) * (maxB - minB);
    const int countfromreferencespheres = 4;
    const int total = count + countfromreferencespheres;
    struct Sphere *world = malloc(total * sizeof(struct Sphere));

    world[0] = groundSphere;
    world[1] = referenceSphere1;
    world[2] = referenceSphere2;
    world[3] = referenceSphere3;
    int index = 4;

    for (int a = minA; a < maxA; a++)
    {
        for (int b = minB; b < maxB; b++)
        {
            float choose_mat = randomFloat(0, 1);
            float center[3];
            center[0] = a + 0.9 * randomFloat(0, 1);
            center[1] = 0.2;
            center[2] = b + 0.9 * randomFloat(0, 1);
            float intermediatestep1[3];
            float test[3] = {4, 0.2, 0};
            subtractVectors(center, test, intermediatestep1);
            float distance = vectorLength(intermediatestep1);

            if (distance > 0.9)
            {
                struct Sphere sphere;
                if (choose_mat < 0.8)
                {
                    sphere.center[0] = center[0];
                    sphere.center[1] = center[1];
                    sphere.center[2] = center[2];
                    sphere.radius = 0.2;
                    sphere.material.type = LAMBERTIAN;
                    randomVector3(0, 1, sphere.material.albedo);
                    world[index] = sphere;
                    index++;
                }
                else if (choose_mat < 0.95)
                {
                    sphere.center[0] = center[0];
                    sphere.center[1] = center[1];
                    sphere.center[2] = center[2];
                    sphere.radius = 0.2;
                    sphere.material.type = METAL;
                    randomVector3(0.5, 1, sphere.material.albedo);
                    sphere.material.fuzz = randomFloat(0, 0.5);
                    world[index] = sphere;
                    index++;
                }
                else
                {
                    sphere.center[0] = center[0];
                    sphere.center[1] = center[1];
                    sphere.center[2] = center[2];
                    sphere.radius = 0.2;
                    sphere.material.type = GLASS;
                    sphere.material.albedo[0] = 1;
                    sphere.material.albedo[1] = 1;
                    sphere.material.albedo[2] = 1;
                    sphere.material.ir = 1.5;
                    world[index] = sphere;
                    index++;
                }
            }
        }
    }

    FILE *f = fopen("test.ppm", "w");
    fprintf(f, "P3\n%i %i\n%i\n", imageWidth, imageHeight, maxVal);
    for (int j = 0; j < imageHeight; j++)
    {
        printf("Scanlines remaining: %i\n", imageHeight - j);
        for (int i = 0; i < imageWidth; i++)
        {
            float color[3] = {0, 0, 0};
            struct Ray r = getRay(i, j, pixel00Location, pixelDeltaU, pixelDeltaV, cameraCenter);
            float raycolor[3];
            for (int k = 0; k < samplesPerPixel; k++)
            {
                rayColor(r, maxDepth, raycolor, world, total);
                addVectors(color, raycolor, color);
            }
            writeColor(f, color, samplesPerPixel);
        }
    }
    return 0;
}
