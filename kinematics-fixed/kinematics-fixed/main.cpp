#define _CRT_SECURE_NO_WARNINGS 
#include "kinematics_fixed.h"
#include "m_mcpy.h"
#include "utils.h"
#include <stdio.h>
#include "dh_hex_fixed.h"
#include <stdio.h>

void main()
{
//	FILE* fp;
//	fp = fopen("log.txt", "W");

	dynamic_hex_t dh_f;
	setup_dynamic_hex(&dh_f);

	for(int leg = 0; leg < 6; leg++)
	{
		for (int joint = 0; joint < 3; joint++)
		{
			joint32_t* j = &(dh_f.p_joint[leg][joint]);
			j->q = 0;
			j->cos_q = cos_lookup(j->q, j->n_r);
			j->sin_q = sin_lookup(j->q, j->n_r);
		}
		forward_kinematics_64(&dh_f.hb_0[leg], dh_f.p_joint[leg]);
	}
	printf("h0_1 = \r\n");
	print_mat4_32b(dh_f.p_joint[0]->h_im1_i);
	printf("\nh1_2 = \n");
	print_mat4_32b(dh_f.p_joint[0]->child->h_im1_i);
	printf("\nh2_3 = \n");
	print_mat4_32b(dh_f.p_joint[0]->child->child->h_im1_i);
	printf("\n\n");

	printf("hb_1 = \r\n");
	print_mat4_32b(dh_f.p_joint[0]->hb_i);
	printf("\nhb_2 = \n");
	print_mat4_32b(dh_f.p_joint[0]->child->hb_i);
	printf("\nhb_3 = \n");
	print_mat4_32b(dh_f.p_joint[0]->child->child->hb_i);
	printf("\n\n");
}

