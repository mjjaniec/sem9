#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

//const double SUN_MASS   = 1.9855e+30;
//const double LIGHT_YEAR = 9.4607e+15;
//const double GRAVITY    = 6.6738e-11;

const double SUN_MASS   = 10;
const double LIGHT_YEAR = 1;
const double GRAVITY    = 1;
const double PI_CONST   = 3.14159265;

double math_square(double x) {
    return x * x;
}

double math_qube(double x) {
    return x * x * x;
}

double math_abs(double x) {
    return x < 0 ? -x : x;
}

typedef struct __V3 {
    double x, y, z;
}V3;

bool V3_isZero(V3 v) {
    return math_abs(v.x) + math_abs(v.y) + math_abs(v.z) < 1e-6;
}

V3 V3_create(double x, double y, double z) {
    V3 result;
    result.x = x;
    result.y = y;
    result.z = z;
    return result;
}

V3 V3_sub(V3 a, V3 b) {
    return V3_create(a.x - b.x, a.y - b.y, a.z - b.z);
}

V3 V3_add(V3 a, V3 b) {
    return V3_create(a.x + b.x, a.y + b.y, a.z + b.z);
}

V3 V3_add_assign(V3 a, V3 b) {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

V3 V3_mul_sv(double s, V3 v) {
    return V3_create(s * v.x, s * v.y, s * v.z);
}

bool V3_equals(V3 a, V3 b) {
    return V3_isZero(V3_sub(a, b));
}

double V3_norm(V3 v){
    return sqrt(math_square(v.x) + math_square(v.y) + math_square(v.z));
}

void V3_print(V3 v, bool nl) {
    printf("[%+5f, %+5f, %+5f]", (float)v.x, (float)v.y, (float)v.z);
    if (nl) {
        printf("\n");
    } else {
        printf("   ");
    }
}

typedef struct __Start {
    V3 position;
//    V3 speed;
    V3 force;
    double mass;
}Star;


double rnd_0_1() {
    return (double)rand()/(double)(RAND_MAX);
}

double rnd_exp(double mean) {
    return -log(rnd_0_1()) * mean;
}



void Star_init(Star* self) {
    self->mass = rnd_exp(1) * SUN_MASS;

    double radius = 20 * LIGHT_YEAR;
    double phi = rnd_0_1()*2*PI_CONST;
    double theta = (rnd_0_1() - 0.5) * PI_CONST / 10;
    self->position = V3_create(radius * sin(phi) * cos(theta),
                     radius * cos(phi) * cos(theta),
                     radius * sin(theta));

//    phi += PI_CONST / 2;
//    theta = -theta / (2 + rnd_0_1());
//    radius /= 20 * (2 + rnd_0_1());
//    self->speed = V3_create(radius * sin(phi)*cos(theta),
//                           radius * cos(phi)*cos(theta),
//                           radius * cos(theta) / 10.0);

    self->force = V3_create(0.0, 0.0, 0.0);
}

void Star_print(Star* self) {
    printf("Star: {\n    mass: %+6.3e\n    pos:   [%+6.3e, %+6.3e, %+6.3e]\n    force: [%+6.3e, %+6.3e, %+6.3e]\n}\n",
           (float)self->mass,
           (float)self->position.x, (float)self->position.y, (float)self->position.z,
  //         (float)self->speed.x, (float)self->speed.y, (float)self->speed.z,
           (float)self->force.x, (float)self->force.y, (float)self->force.z);
}
