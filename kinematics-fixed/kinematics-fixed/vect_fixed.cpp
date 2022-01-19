/*
 * vect.c
 *
 *  Created on: Nov 21, 2021
 *      Author: Ocanath Robotman
 * 
 * 
 * 
*	Square Root methods taken DIRECTLY FROM: 
*	https://github.com/chmike/fpsqrt/blob/master/fpsqrt.c
*	Credit for square root algorithms goes to Christophe Meessen
* 
 */
#include "vect_fixed.h"

#define ONE_12B 4096

 /*
  * Same as normalt mat4 multiplication with pass by pointer,
  * EXCEPT that it is assummed to be a homogeneous transformation matrix where
  * the bottom rows of all matrices are equal to 0,0,0,1
 */
void ht32_mult_pbr(mat4_32b_t* m1, mat4_32b_t* m2, mat4_32b_t* ret)
{
	int r, c, i;
	int32_t tmp;
	/*Do the rotation matrix part FIRST*/
	for (r = 0; r < 3; r++)
	{
		for (c = 0; c < 3; c++)
		{
			tmp = 0;
			for (i = 0; i < 4; i++)
			{
				tmp += (m1->m[r][i] * m2->m[i][c]);
			}
			ret->m[r][c] = tmp >> 12;	//each term has 12bit binary decimal point. Remove the additional terms accumulated AFTER the sum.
			//the 3x3 portion is safe to do multiple accumulations because each term is 12bit limited (even ignoring constraints imposed by the
			//conditions required to be a homogeneous transformation matrix, such as column vector orthogonality, unit vector across row and column vectors)
		}
	}
	/*Do the translation part LAST. This part must be handled differently due to relative scaling*/

	for (int r = 0; r < 4; r++)
	{
		tmp = 0;
		for (int i = 0; i < 4; i++)
			tmp += m1->m[r][i] * m2->m[i][3];
		ret->m[r][3] = tmp >> 12;	//remove the order 12 component
	}

	//	ret->m[3][0] = 0.f;
	//	ret->m[3][1] = 0.f;
	//	ret->m[3][2] = 0.f;
	//	ret->m[3][3] = 1.f;	//commented because this function should be called on matrices where that relationship is assumed to be true
}

/*
 * Same as normalt mat4 multiplication with pass by pointer,
 * EXCEPT that it is assummed to be a homogeneous transformation matrix where
 * the bottom rows of all matrices are equal to 0,0,0,1
*/
void ht32_mult64_pbr(mat4_32b_t* m1, mat4_32b_t* m2, mat4_32b_t* ret, int n)
{
	//if(n > 29)
	//	return
	int r, c, i;
	int64_t tmp;
	/*Do the rotation matrix part FIRST*/
	for (r = 0; r < 4; r++)
	{
		for (c = 0; c < 4; c++)
		{
			tmp = 0;
			for (i = 0; i < 4; i++)
			{
				tmp += (((int64_t)m1->m[r][i]) * ((int64_t)m2->m[i][c]));
			}
			ret->m[r][c] = (int32_t)(tmp >> n);	//shifting after is ok ONLY if the rotation order is less than 29 bits
		}
	}
}

/**/
void  h32_v32_mult(mat4_32b_t* m, vect3_32b_t* vin, vect3_32b_t* vout, int rshift)
{
	int64_t one = ((int64_t)1) << rshift;
	for (int r = 0; r < 3; r++)
	{
		int64_t tmp = 0;
		for (int c = 0; c < 3; c++)
		{
			tmp += ( ((int64_t)m->m[r][c]) * ((int64_t)vin->v[c]) );
		}
		tmp += (int64_t)m->m[r][3]*one;	
		vout->v[r] = (int32_t)(tmp >> rshift);
	}
}

/*
	Returns vector cross product between 3 vectors A and B. Faster pass by pointer version
*/
void cross32_pbr(vect3_32b_t* v_a, vect3_32b_t* v_b, vect3_32b_t* ret, int n)
{
	ret->v[0] = (-v_a->v[2] * v_b->v[1] + v_a->v[1] * v_b->v[2]) >> n;
	ret->v[1] = (v_a->v[2] * v_b->v[0] - v_a->v[0] * v_b->v[2]) >> n;
	ret->v[2] = (-v_a->v[1] * v_b->v[0] + v_a->v[0] * v_b->v[1]) >> n;	//need to remove some of the scaling factor used to express these in fixed point
}

/*
	Same as cross32, except that a 64bit buffer is used, allowing full scale
	resolution for input vectors.
*/
void cross64_pbr(vect3_32b_t* v_a, vect3_32b_t* v_b, vect3_32b_t* ret, int rshift)
{
	//ret->v[0] = (-v_a->v[2] * v_b->v[1] + v_a->v[1] * v_b->v[2]) >> n;
	//ret->v[1] = (v_a->v[2] * v_b->v[0] - v_a->v[0] * v_b->v[2]) >> n;
	//ret->v[2] = (-v_a->v[1] * v_b->v[0] + v_a->v[0] * v_b->v[1]) >> n;	//need to remove some of the scaling factor used to express these in fixed point
	int64_t va0 = (int64_t)v_a->v[0];
	int64_t va1 = (int64_t)v_a->v[1];
	int64_t va2 = (int64_t)v_a->v[2];

	int64_t vb0 = (int64_t)v_b->v[0];
	int64_t vb1 = (int64_t)v_b->v[1];
	int64_t vb2 = (int64_t)v_b->v[2];

	ret->v[0] = (int32_t)((-va2 * vb1 + va1 * vb2) >> rshift);
	ret->v[1] = (int32_t)((va2 * vb0 - va0 * vb2) >> rshift);
	ret->v[2] = (int32_t)((-va1 * vb0 + va0 * vb1) >> rshift);
}

/*Helper function to convert array of floats to array of int32 with proper radix representation*/
void float_to_int32(float* in, int32_t* out, int num_elements, int radix)
{
	float c = (float)(1 << radix);	//constant to multiply up input in order to get proper radix
	for (int i = 0; i < num_elements; i++)
		out[i] = (int32_t)(in[i] * c);
}

int32_t vect64_mag(vect3_32b_t* v)
{
	return 0;
}


// sqrt_i32 computes the squrare root of a 32bit integer and returns
// a 32bit integer value. It requires that v is positive.
//longer compute time for bigger inputs v
int32_t sqrt_i32(int32_t v) 
{
	uint32_t b = 1 << 30, q = 0, r = v;
	while (b > r)
		b >>= 2;
	while (b > 0) 
	{
		uint32_t t = q + b;
		q >>= 1;
		if (r >= t) 
		{
			r -= t;
			q += b;
		}
		b >>= 2;
	}
	return q;
}

// sqrt_i64 computes the squrare root of a 64bit integer and returns
// a 64bit integer value. It requires that v is positive.
int64_t sqrt_i64(int64_t v) 
{
	uint64_t b = ((uint64_t)1) << 62, q = 0, r = v;
	while (b > r)
		b >>= 2;
	while (b > 0) 
	{
		uint64_t t = q + b;
		q >>= 1;
		if (r >= t)
		{
			r -= t;
			q += b;
		}
		b >>= 2;
	}
	return q;
}


int64_t vect32_mag64(vect3_32b_t* vin, int vin_radix)
{
	int64_t varr[3];	//64bit representation of vin
	int64_t v_dot_v = 0;
	for (int i = 0; i < 3; i++)
	{
		varr[i] = (int64_t)vin->v[i];
		v_dot_v += (varr[i] * varr[i]) >> vin_radix;
	}
	//vdotv is guaranteed positive.
	return  sqrt_i64(v_dot_v << vin_radix);	//left shift by the radix because sqrt(V*2^32) = 2^16*sqrt(v). Pre-left shifting preserves the radix of the input without loss of resolution
}

/*
* Expensive. 64 bit buffered fixed point dot product and 
* inverse square root. Tested, works. No guarding
* for overflow, so be careful...
* 
* 
* Note. Repeating the functionality of vect32_mag64 here wihtout
* an actual function call because both sections require the use of
* a 64bit buffered copy of vin. This is a slight optimization
* at the expense of clutter/code re-use
*/
void normalize_vect64(vect3_32b_t* vin, int vin_radix)
{
	//obtain the magnitude
	int64_t varr[3];	//64bit representation of vin
	int64_t v_dot_v = 0;	//note: vdotv is guaranteed positive.
	for (int i = 0; i < 3; i++)
	{
		varr[i] = (int64_t)vin->v[i];
		v_dot_v += (varr[i] * varr[i]) >> vin_radix;	
	}
	int64_t mag = sqrt_i64(v_dot_v << vin_radix);	//left shift by the radix because sqrt(V*2^32) = 2^16*sqrt(v). Pre-left shifting preserves the radix of the input without loss of resolution

	//scale the input vector
	int64_t one = ((int64_t)1 << vin_radix);	//t y p e s
	for(int i = 0; i < 3; i++)
		vin->v[i] = (int32_t)((varr[i] * one) / mag);	//apply 1/sq. Multiply by radix to preserve sign, and premultiply by one to prevent truncation of lower bits. mag is same radix as vin
}

/*Generic 64bit buffer dot product calculator for vectors of arbitrary size*/
int32_t dot64_pbr(int32_t* va, int32_t* vb, int size, int rshift)
{
	int64_t tmp = 0;
	for (int i = 0; i < size; i++)
	{
		int64_t va64 = (int64_t)va[i];
		int64_t vb64 = (int64_t)vb[i];
		tmp += ((va64 * vb64) >> rshift);	//could wait until the end to rshift, but there could be overflow if the sum of the multiples exceeds 64 bits
	}
	return (int32_t)tmp;
}

void h32_origin_pbr(vect3_32b_t * v, mat4_32b_t* m)
{
	for (int r = 2; r >= 0; r--)
		v->v[r] = m->m[r][3];
}

vect3_32b_t h32_origin(mat4_32b_t* m)
{
	vect3_32b_t ret;
	for (int r = 0; r < 3; r++)
		ret.v[r] = m->m[r][3];
	return ret;
}

/*Loads rotation about coordinate. 0 = identity*/
mat4_32b_t Hz_nb(int32_t angle, int n)
{
	int32_t cth = cos_lookup(angle, n);
	int32_t sth = sin_lookup(angle, n);
	int32_t one_n = (1 << n);
	mat4_32b_t r;
	r.m[0][0] = cth;		r.m[0][1] = -sth;		r.m[0][2] = 0;			r.m[0][3] = 0;
	r.m[1][0] = sth;		r.m[1][1] = cth;		r.m[1][2] = 0;			r.m[1][3] = 0;
	r.m[2][0] = 0;			r.m[2][1] = 0;			r.m[2][2] = one_n;		r.m[2][3] = 0;
	r.m[3][0] = 0;			r.m[3][1] = 0;			r.m[3][2] = 0;			r.m[3][3] = one_n;
	return r;
}

/*returns homogeneous transform mat4_t matrix which is the rotation 'analge' around the x axis */
mat4_32b_t Hx_nb(int32_t angle, int n)
{
	int32_t cth = cos_lookup(angle, n);
	int32_t sth = sin_lookup(angle, n);
	int32_t one_n = (1 << n);
	mat4_32b_t ret;
	ret.m[0][0] = one_n;	ret.m[0][1] = 0;	ret.m[0][2] = 0;		ret.m[0][3] = 0;
	ret.m[1][0] = 0;		ret.m[1][1] = cth;	ret.m[1][2] = -sth;		ret.m[1][3] = 0;
	ret.m[2][0] = 0;		ret.m[2][1] = sth;	ret.m[2][2] = cth;		ret.m[2][3] = 0;
	ret.m[3][0] = 0;		ret.m[3][1] = 0;	ret.m[3][2] = 0;		ret.m[3][3] = one_n;
	return ret;
}

/*Returns rotation about coordinate. 0 = identity*/
mat4_32b_t Hy_nb(int32_t angle, int n)
{
	int32_t cth = cos_lookup(angle, n);
	int32_t sth = sin_lookup(angle, n);
	int32_t one_n = (1 << n);
	mat4_32b_t ret;
	ret.m[0][0] = cth;		ret.m[0][1] = 0;		ret.m[0][2] = sth;		ret.m[0][3] = 0;
	ret.m[1][0] = 0;		ret.m[1][1] = one_n;	ret.m[1][2] = 0;		ret.m[1][3] = 0;
	ret.m[2][0] = -sth;		ret.m[2][1] = 0;		ret.m[2][2] = cth;		ret.m[2][3] = 0;
	ret.m[3][0] = 0;		ret.m[3][1] = 0;		ret.m[3][2] = 0;		ret.m[3][3] = one_n;
	return ret;
}

