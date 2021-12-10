#include "utils.h"
#include <stdio.h>


void print_mat4_32b(mat4_32b_t m)
{
	for (int r = 0; r < 4; r++)
	{

		for (int c = 0; c < 4; c++)
		{
			float v;
			if (c == 3 && r < 3)
				v = ((float)m.m[r][c]) / ((float)(1 << 16));
			else
				v = ((float)m.m[r][c]) / ((float)(1 << 21));
			printf("%f ", v);
		}
		printf("\n");
	}
}
