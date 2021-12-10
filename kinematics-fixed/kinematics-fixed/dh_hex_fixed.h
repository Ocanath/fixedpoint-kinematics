/*
 * hexapod_params.h
 *
 *  Created on: Nov 21, 2021
 *      Author: Ocanath Robotman
 */

#ifndef INC_HEXAPOD_PARAMS_H_
#define INC_HEXAPOD_PARAMS_H_
#include "kinematics_fixed.h"

#define NUM_LEGS			6
#define NUM_DOFS_PER_LEG	3
#define NUM_JOINTS_TOTAL	NUM_LEGS*NUM_DOFS_PER_LEG

typedef struct dynamic_hex_t
{
	mat4_32b_t hb_0[NUM_LEGS];
	joint32_t * p_joint[NUM_LEGS];
}dynamic_hex_t;

extern dynamic_hex_t gl_hex;

void setup_dynamic_hex(dynamic_hex_t* robot);


#endif /* INC_HEXAPOD_PARAMS_H_ */
