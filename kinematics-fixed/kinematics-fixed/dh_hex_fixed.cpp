#include "dh_hex_fixed.h"
#include "m_mcpy.h"

joint32_t chain[NUM_DOFS_PER_LEG * NUM_LEGS];	//the memory used to describe the robot. the robot struct contains pointers to this memory

static const mat4_32b_t hb_0 = {
	{    //hb_0
		{-ONE_ROT, 0, 0, 7194922},
		{0, ONE_ROT, 0, 0},
		{0, 0, -ONE_ROT, 0},
		{0, 0, 0, ONE_ROT}
	}
};
/*Optional: compute these with urdf-style rotations and translations
* or DH entries. These were pre-computed in matlab for this example.
*/
static const mat4_32b_t links_def[NUM_LEGS + 1] = {
	{
		{    //link 1
			{ONE_ROT, 0, 0, -3486515 },
			{0, 0, -ONE_ROT, 0},
			{0, ONE_ROT, 0, 4303094 },
			{0, 0, 0, ONE_ROT}
		}
	},
	{
		{    //link 2
			{ONE_ROT, 0, 0, -6584141},
			{0, -ONE_ROT, 0, 0},
			{0, 0, -ONE_ROT, 1900544},
			{0, 0, 0, ONE_ROT}
		}
	},
	{
		{    //link 3
			{ONE_ROT, 0, 0, -12996888},
			{0, ONE_ROT, 0, 0},
			{0, 0, ONE_ROT, 1409024},
			{0, 0, 0, ONE_ROT}
		}
	}
};

void setup_dynamic_hex(dynamic_hex_t* robot)
{
	/*Set up the chain singly linked list. This formats the chain
	 * in ascending order, with each triple corresponding
	 * to a leg. */
	int leg = 0;
	int last_idx = NUM_JOINTS_TOTAL - 1;
	for (int i = 0; i < NUM_JOINTS_TOTAL; i++)
	{
		//lookup the current element once. Also keeps code less busy
		joint32_t* j = &chain[i];	
		
		//Load the frame number of this joint
		int frame = (i % NUM_DOFS_PER_LEG) + 1;
		j->frame = frame;

		//Initialize the labelled orders of rotations, translations, and jacobian
		j->n_r = KINEMATICS_SIN_ORDER;
		j->n_t = KINEMATICS_TRANSLATION_ORDER;
		j->n_si = KINEMATICS_SIN_ORDER;
		
		j->q = 0;
		j->sin_q = sin_lookup(j->q, KINEMATICS_SIN_ORDER);
		j->cos_q = cos_lookup(j->q, KINEMATICS_SIN_ORDER);

		//load the 'variable' matrices with identity, and the link frames from the initializer list
		m_mcpy(&j->hb_i, (void*)gl_identity_matrix_32b, sizeof(mat4_32b));
		m_mcpy(&j->h_im1_i, (void*)gl_identity_matrix_32b, sizeof(mat4_32b));
		m_mcpy(&j->h_link, (void*)(&links_def[frame-1]), sizeof(mat4_32b_t));

		//Establish the singly linked list
		if(frame < 3)
			j->child = &chain[i + 1];
		else
			j->child = NULL;	//terminate the SLL with a null pointer. 

		if (frame == 1)
			robot->p_joint[leg++] = &chain[i];
	}

	/*Initialize the hb_0 frame definitions*/
	const int32_t angle_12 = TWO_PI_12B / 6;
	for (int leg = 0; leg < NUM_LEGS; leg++)
	{
		int n = KINEMATICS_SIN_ORDER;

		mat4_32b_t* rh32_b_0 = &(robot->hb_0[leg]);
		
		mat4_32b_t z32_rot = Hz_nb(leg * angle_12, KINEMATICS_SIN_ORDER);	//get the z rotationz

		ht32_mult64_pbr(&z32_rot, (mat4_32b_t*)(&hb_0), rh32_b_0, KINEMATICS_SIN_ORDER);
	}
}