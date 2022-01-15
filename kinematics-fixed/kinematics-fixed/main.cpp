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

void print_vect_mm(const char* prefix, vect3_32b_t* v, int radix)
{
	float div = (float)(1 << radix);
	float res[3];
	for (int i = 0; i < 3; i++)
		res[i] = (float)v->v[i] / div;
	printf("%s:[%f,%f,%f]\r\n", prefix, res[0], res[1], res[2]);
}


int gradient_descent_ik(mat4_32b_t * hb_0, joint32_t* start, joint32_t* end, vect3_32b_t* o_targ_b, vect3_32b_t* o_anchor_b)
{
	int solved = 0;
	print_vect_mm("TARG", o_targ_b, start->n_t);
	int cycles = 0;
	while (solved == 0)
	{
		//do forward kinematics
		forward_kinematics_64(hb_0, start);
		h32_origin_pbr(o_anchor_b, &end->hb_i);
		calc_J_32b_point(hb_0, start, o_anchor_b);

		//show progress
		printf("targ: [%d,%d,%d], anchor pos: [%d,%d,%d]\r\n", o_targ_b->v[0], o_targ_b->v[1], o_targ_b->v[2], o_anchor_b->v[0], o_anchor_b->v[1], o_anchor_b->v[2]);
		//print_vect_mm("anchor", &o_anchor_b, j->n_t);
		//printf("q: [%f, %f, %f]\r\n", (float)j[0].q / PI_12B_BY_180, (float)j[1].q / PI_12B_BY_180, (float)j[2].q / PI_12B_BY_180);

		//get vector pointing from the anchor point on the robot to the target. call it 'f'
		vect3_32b_t f;
		for (int i = 0; i < 3; i++)
			f.v[i] = (o_targ_b->v[i] - o_anchor_b->v[i]) / 500;

		//get the static torque produced by the force vector. Radix should be same as established in 'f' if j->si
		int tau_rshift = start->n_t;
		calc_j_taulist(start, &f, tau_rshift);	//removing an n_si (from f) yields tau in radix 16
		int tau_radix = (start->n_si + start->n_t) - tau_rshift;

		//iterate through joints. could be sll traversal
		int32_t one = 1 << start->n_r;
		vect3_32b_t z = { 0, 0, one };
		solved = 1;
		joint32_t* j = start;
		while (j != NULL)
		{
			vect3_32b_t vq = { j->cos_q, j->sin_q, 0 };

			vect3_32b_t tangent;
			cross64_pbr(&z, &vq, &tangent, j->n_r);
			vect3_32b_t vq_new;

			int64_t tau_i_64 = (int64_t)j->tau_static;
			for (int r = 0; r < 3; r++)
			{
				int64_t tmp = (((int64_t)tangent.v[r]) * tau_i_64) >> tau_radix;
				tmp /= 150;
				vq_new.v[r] = (int32_t)tmp + vq.v[r];

				if (tmp != 0)
					solved = 0;
			}

			normalize_vect64(&vq_new, j->n_si);

			j->cos_q = vq_new.v[0];
			j->sin_q = vq_new.v[1];

			j = j->child;
		}
		cycles++;
	}
	int rshift = (start->n_r - 12);
	joint32_t* j = start;
	while (j != NULL)
	{
		j->q = atan2_fixed(j->sin_q >> rshift, j->cos_q >> rshift);
		j = j->child;
	}
	return cycles;
}

int main(void)
{
	{
		int64_t v = 84 * (1 << KINEMATICS_TRANSLATION_ORDER);
		int64_t vnew = sqrt_i64(v);
		printf("in: %d, out: %d\r\n", (int32_t)v, (int32_t)vnew);
	}
	{
		float f[3] = { 4.7, 9.21, -11.15 };
		vect3_32b_t in;
		int n = KINEMATICS_SIN_ORDER;
		float scf = (float)(1 << n);
		for (int i = 0; i < 3; i++)
			in.v[i] = (int32_t)(f[i] * scf);
		print_vect_mm("input: ", &in, n);
		normalize_vect64(&in, n);
		print_vect_mm("normalized: ", &in, n);
	}

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
		o_footip_3.v[i] = (int32_t)(o_footip_3_f[i] * start->n_t);	//establish a fixed-point representation of the foot tip position in mm*translationradix, in coordinate frame 3

	/*
	* Note: for convenience, the singly linked list elements are ALSO in an array.
	* However, the indexing is now broken. j[2] is frame 3.
	*/
	start[0].q = fdeg_to_12b(80.f);	//q1
	start[1].q = fdeg_to_12b(-30.f);	//q2
	start[2].q = fdeg_to_12b(-15.f);	//q3
	int32_t qtarg[3];
	for (int i = 0; i < 3; i++)
		qtarg[i] = start[i].q;
	load_qsin(start);
	forward_kinematics_64(m, start);
	
	//vect3_32b_t otarg = h32_origin_pbr(&otarg, j[2].hb_i);
	vect3_32b_t otarg;
	for (int i = 0; i < 3; i++)
		otarg.v[i] = start[2].hb_i.m[i][3];	//otarg is a point in the base frame
	vect3_32b_t o_anchor_b;	//represents the anchor point for the gradient descent force vector


	/*Re-initialize the joint positions to 0*/
	start[0].q = fdeg_to_12b(0.f);
	start[1].q = fdeg_to_12b(-100.f);
	start[2].q = fdeg_to_12b(-100.f);
	load_qsin(start);

	int cycles = gradient_descent_ik(m, start, end, &otarg, &o_anchor_b);

	float div = (float)(1 << KINEMATICS_TRANSLATION_ORDER);
	float res[3];
	for (int i = 0; i < 3; i++)
		res[i] = (float)(otarg.v[i] - o_anchor_b.v[i])/div;
	
	printf("Final Error: [%f,%f,%f]\r\n", res[0], res[1], res[2]);
	printf("Computed in %d cycles\r\n", cycles);
	div = 71.4887f;
	printf("Q: [%f,%f,%f]\r\n", (float)start[0].q / div, (float)start[1].q / div, (float)start[2].q / div);
}