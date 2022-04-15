#include "print_state.h"

void CAL_CONV gfexfo_(int* iwa, float* rwa, double* dwa, int* ipdms,
	int* iinfo, float* rinfo,
	int* npcur, float* curcor, float* curvel,
	int* kxfo, float* rxfo, float* rhxfo,
	int* ixfo, int* iextf, int* icoord,
	int* nint, int* nrea, int* nsto, int* nstr,
	char chext[][MCHEXT], float state[][MSTATE],
	float gamma[][MGAMMA], float vellot[][MVELLOT],
	float* stor, int* ierr)
{
	int step_nr = iinfo[5];
	int substep_nr = iinfo[7];
	int iteration_nr = iinfo[11];

	//if (step_nr == iinfo[6] && substep_nr == iinfo[8] && iteration_nr == 0)
	{
		int body_number=iinfo[2]-1;
		int number_of_bodies=iinfo[4];
		
		//for (body_number=0;body_number<number_of_bodies;body_number++)
		{
			//printf("State array for body %d out of %d:\n",body_number,number_of_bodies);
			for (int i = 1; i<11; i=i+2)
			{
				printf("%9.1E",state[body_number][i]);
			}
			printf("\n");
		}
	}

	for(int i = 0; i < *nsto; i++)
	{
		rhxfo[i] = 0.0; // Store zero forces in intermediate array
	}
	for(int i = 0; i < 9; i++)
	{
		stor[i] = 0.; // Store zero forces in result array
		stor[0] = 5000.; // [kN] surge force
	}

	*ierr = 0;
}

