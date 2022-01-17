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

#define PI_12B_BY_180	71.488686f

int32_t fdeg_to_12b(float f)
{
	return (int32_t)(f * PI_12B_BY_180);	//4096 * pi / 180 = 71.488686
}

static const float o_footip_3_f[3] = { -15.31409f, -9.55025f, 0.f };

void print_vect_mm(const char* prefix, vect3_32b_t* v, int radix, const char * suffix)
{
	float div = (float)(1 << radix);
	float res[3];
	for (int i = 0; i < 3; i++)
		res[i] = (float)v->v[i] / div;
	printf("%s:[%f,%f,%f]%s", prefix, res[0], res[1], res[2], suffix);
}

int main(void)
{
	std::ofstream fp_position;
	fp_position.open("ik_efpos.csv");
	std::ofstream fp_q;
	fp_q.open("ik_q.csv");

	dynamic_hex_t dh_f;
	setup_dynamic_hex(&dh_f);	//set up a dynamic hexapod robot structure
	
	int leg = 0;
	//joint32_t * j = dh_f.p_joint[leg];
	joint32_t * start = dh_f.p_joint[leg];
	joint32_t* end = &start[2];
	mat4_32b_t* m = &dh_f.hb_0[leg];
	vect3_32b_t o_footip_3;
	for (int i = 0; i < 3; i++)
		o_footip_3.v[i] = (int32_t)(o_footip_3_f[i] * (float)(1<<start->n_t) );	//establish a fixed-point representation of the foot tip position in mm*translationradix, in coordinate frame 3

	/*
	* Note: for convenience, the singly linked list elements are ALSO in an array.
	* However, the indexing is now broken. j[2] is frame 3.
	*/
	start[0].q = fdeg_to_12b(15.f);	//q1
	start[1].q = fdeg_to_12b(-15.f);	//q2
	start[2].q = fdeg_to_12b(90.f);	//q3
	int32_t qtarg[3];
	for (int i = 0; i < 3; i++)
		qtarg[i] = start[i].q;
	load_qsin(start);
	forward_kinematics_64(m, start);
	
	//vect3_32b_t otarg = h32_origin_pbr(&otarg, j[2].hb_i);
	vect3_32b_t ostart;
	h32_v32_mult(&start[2].hb_i, &o_footip_3, &ostart, start->n_r);
	
	/*Re-initialize the joint positions to 0*/
	start[0].q = fdeg_to_12b(0.f);
	start[1].q = fdeg_to_12b(0.f);
	start[2].q = fdeg_to_12b(0.f);
	load_qsin(start);

	float one_t = (float)(1 << start->n_t);
	for (float x = -30.f; x < 30.f; x += 3.f)
	{
		for (float y = -30.f; y < 30.f; y += 3.f)
		{
			for (float z = -30.f; z < 30.f; z += 3.f)
			{
				vect3_32b_t otarg;
				float arr[] = { x,y,z };
				for (int i = 0; i < 3; i++)
					otarg.v[i] = (int32_t)(arr[i] * one_t)+ostart.v[i];

				print_vect_mm("targ:", &otarg, start->n_t, " ");
				vect3_32b_t o_anchor_b;	//represents the anchor point for the gradient descent force vector
				int cycles = gradient_descent_ik(m, start, end, &o_footip_3, &otarg, &o_anchor_b, 7000);
				print_vect_mm("anchor:", &o_anchor_b, start->n_t, " ");
				printf("cycles: %d, ", cycles);
				float div = 71.4887f;
				printf("Q: [%f,%f,%f]\r\n", (float)start[0].q / div, (float)start[1].q / div, (float)start[2].q / div);
			}
		}
	}
}