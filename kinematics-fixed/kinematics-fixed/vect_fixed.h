#ifndef VECT_FIXED_H
#define VECT_FIXED_H
#include "trig_fixed.h"

typedef int32_t mat4_32b[4][4];
typedef struct mat4_32b_t
{
	mat4_32b m;
}mat4_32b_t;

typedef int32_t vect3_32b[3];
typedef struct vect3_32b_t
{
	vect3_32b v;
}vect3_32b_t;
typedef int32_t vect6_32b[6];
typedef struct vect6_32b_t
{
	vect6_32b v;
}vect6_32b_t;

typedef int32_t quat_32b[4];
typedef struct quat_32b_t
{
	quat_32b v;
}quat_32b_t;	//quat is a 4 vector

mat4_32b_t Hy_nb(int32_t angle, int n);
mat4_32b_t Hx_nb(int32_t angle, int n);
mat4_32b_t Hz_nb(int32_t angle, int n);
void h32_origin_pbr(vect3_32b_t* v, mat4_32b_t* m);

void ht32_mult_pbr(mat4_32b_t* m1, mat4_32b_t* m2, mat4_32b_t* ret);
void cross32_pbr(vect3_32b_t* v_a, vect3_32b_t* v_b, vect3_32b_t* ret, int n);
void cross64_pbr(vect3_32b_t* v_a, vect3_32b_t* v_b, vect3_32b_t* ret, int n);
void ht32_mult64_pbr(mat4_32b_t* m1, mat4_32b_t* m2, mat4_32b_t* ret, int n);

#endif // !VECT_FIXED_H
