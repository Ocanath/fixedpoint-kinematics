#define _CRT_SECURE_NO_WARNINGS 
#include "kinematics_fixed.h"
#include "m_mcpy.h"
#include "utils.h"
#include <stdio.h>
#include "dh_hex_fixed.h"

#include <iostream>
#include <fstream>

void create_workspace_file(void)
{
	std::ofstream fp_position;
	fp_position.open("sweep_fixed.csv");
	std::ofstream fp_q;
	fp_q.open("sweep_q.csv");

	dynamic_hex_t dh_f;
	setup_dynamic_hex(&dh_f);
	
	int r = 10;
	for (int i = 0; i <= 0; i++)
	{
		for (int ji = -r; ji <= r; ji++)
		{
			for(int k = -r; k <= r; k++)
			{
				int leg = 0;

				joint32_t* j = dh_f.p_joint[leg];
				j[0].q = (PI_12B * i) / r;
				j[1].q = (PI_12B * ji) / r;
				j[2].q = (PI_12B * k) / r;

				while (j != NULL)
				{
					j->cos_q = cos_lookup(j->q, j->n_r);
					j->sin_q = sin_lookup(j->q, j->n_r);
					j = j->child;
				}

				j = dh_f.p_joint[leg];
				if (j == NULL)
					break;
				mat4_32b_t* m = &dh_f.hb_0[leg];
				forward_kinematics_64(m, j);

				fp_position << dh_f.p_joint[0][2].hb_i.m[0][3] << ", ";
				fp_position << dh_f.p_joint[0][2].hb_i.m[1][3] << ", ";
				fp_position << dh_f.p_joint[0][2].hb_i.m[2][3] << "\n";


				fp_q << j[0].q << ", ";
				fp_q << j[1].q << ", ";
				fp_q << j[2].q << "\n";
			}
		}
	}
	printf("Done\n");

	fp_q.close();
	fp_position.close();
}

int32_t fdeg_to_12b(float f)
{
	return (int32_t)(f * 71.488686f);	//4096 * pi / 180 = 71.488686
}

static const float o_footip_3_f[3] = { -15.31409f, -9.55025f, 0.f };

void print_vect_mm(const char* prefix, vect3_32b_t* v, int radix)
{
	float div = (float)(1 << radix);
	float res[3];
	for (int i = 0; i < 3; i++)
		res[i] = (float)v->v[i] / div;
	printf("%s:[%f,%f,%f]\r\n", prefix, res[0], res[1], res[2]);
}

int main(void)
{
	std::ofstream fp_position;
	fp_position.open("ik_efpos.csv");
	std::ofstream fp_q;
	fp_q.open("ik_q.csv");

	dynamic_hex_t dh_f;
	setup_dynamic_hex(&dh_f);	//set up a dynamic hexapod robot structure
	int32_t tau[3] = { 0 };

	int leg = 0;
	joint32_t * j = dh_f.p_joint[leg];
	mat4_32b_t* m = &dh_f.hb_0[leg];
	vect3_32b_t o_footip_3;
	for (int i = 0; i < 3; i++)
		o_footip_3.v[i] = (int32_t)(o_footip_3_f[i] * j->n_t);	//establish a fixed-point representation of the foot tip position in mm*translationradix, in coordinate frame 3

	/*
	* Note: for convenience, the singly linked list elements are ALSO in an array.
	* However, the indexing is now broken. j[2] is frame 3.
	*/
	j[0].q = fdeg_to_12b(80.f);	//q1
	j[1].q = fdeg_to_12b(-30.f);	//q2
	j[2].q = fdeg_to_12b(-15.f);	//q3
	load_qsin(j);
	forward_kinematics_64(m, j);
	
	//vect3_32b_t otarg = h32_origin_pbr(&otarg, j[2].hb_i);
	vect3_32b_t otarg;
	for (int i = 0; i < 3; i++)
		otarg.v[i] = j[2].hb_i.m[i][3];	//otarg is a point in the base frame
	vect3_32b_t o_anchor_b;	//represents the anchor point for the gradient descent force vector


	/*Re-initialize the joint positions to 0*/
	for (int i = 0; i < 3; i++)
		j[i].q = 0;
	load_qsin(j);

	int solved = 0;
	print_vect_mm("TARG", &otarg, j->n_t);
	while (solved == 0)
	{
		//do forward kinematics
		forward_kinematics_64(m, j);
		h32_origin_pbr(&o_anchor_b, &j[2].hb_i);
		calc_J_32b_point(m, j, &o_anchor_b);

		//get vector pointing from the anchor point on the robot to the target. call it 'f'
		vect3_32b_t f;
		for (int i = 0; i < 3; i++)
			f.v[i] = (otarg.v[i] - o_anchor_b.v[i]);

		//get the static torque produced by the force vector. Radix should be same as established in 'f' if j->si
		calc_j_taulist(j, &f, tau, j->n_si);	//removing an n_t (from f) means tau has radix n_si

		//gradient descent step: increment q to move the anchor in the direction of the target
		solved = 1;
		for (int i = 0; i < 3; i++)
		{
			int32_t step = tau[i] / 850000;
			j[i].q += step;

			if (step != 0)
				solved = 0;	//add solved logic. zero tau means we've stabilized and no changes will occur
		}
		load_qsin(j);

		//show progress
		//printf("anchor pos: [%d,%d,%d]\r\n", o_anchor_b.v[0], o_anchor_b.v[1], o_anchor_b.v[2]);
		print_vect_mm("anchor", &o_anchor_b, j->n_t);
	}
	
	float div = (float)(1 << KINEMATICS_TRANSLATION_ORDER);
	float res[3];
	for (int i = 0; i < 3; i++)
		res[i] = (float)(otarg.v[i] - o_anchor_b.v[i])/div;
	
	printf("Final Error: [%f,%f,%f]\r\n", res[0], res[1], res[2]);
	div = 71.4887f;
	printf("Q: [%f,%f,%f]\r\n", (float)j[0].q / div, (float)j[1].q / div, (float)j[2].q / div);
}