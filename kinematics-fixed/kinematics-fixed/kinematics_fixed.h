#ifndef KINEMATICS_FIXED_H
#define KINEMATICS_FIXED_H
#include "vect_fixed.h"

#define KINEMATICS_SIN_ORDER 			21
#define KINEMATICS_TRANSLATION_ORDER 	16

#define ONE_ROT					(1 << KINEMATICS_SIN_ORDER)

typedef struct joint32_t
{
	int n_t;	//radix, translation part
	int n_r;	//radix, in which sin and cos are expressed (r = rotational part)
	int n_si;	//radix, in which the jacobian column vector Si will be expressed. 

	int frame; 

	mat4_32b_t h_link;
	mat4_32b_t hb_i;
	mat4_32b_t h_im1_i;

	vect6_32b_t Si;

	//TODO: implement and benchmark quaternion maths, see if it is better
	//vect3_32b_t o_link; //origin of the next frame, wrt. our origin.
	//quat_32b_t quat_link;	//orientation of the next frame, expressed in our frame

	int32_t q;
	int32_t sin_q;
	int32_t cos_q;

	int32_t tau_static;

	struct joint32_t* child;
}joint32_t;


extern const mat4_32b gl_identity_matrix_32b;

void forward_kinematics_64(mat4_32b_t* hb_0, joint32_t* f1_joint);
void calc_J_32b_point(mat4_32b_t* hb_0, joint32_t* chain_start, vect3_32b_t* point_b);
void load_qsin(joint32_t* chain_start);
void calc_32b_taulist(joint32_t* chain_start, vect3_32b_t* f, int rshift);
joint32_t* last_joint(joint32_t* chain_start);
int gradient_descent_ik(mat4_32b_t* hb_0, joint32_t* start, joint32_t* end, vect3_32b_t* o_anchor_end, vect3_32b_t* o_targ_b, vect3_32b_t* o_anchor_b, int32_t epsilon_divisor);

#endif // !KINEMATICS_FIXED_H
