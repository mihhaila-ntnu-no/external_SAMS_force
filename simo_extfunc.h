#include <winsock2.h>
#include <ws2tcpip.h>

#include "cal_conv.h"

#pragma comment(lib, "Ws2_32.lib")

#define MSTATE 12
#define MGAMMA 9
#define MCHEXT 120
#define MVELLOT 6

void CAL_CONV gfexfo_(int* iwa, float* rwa, double* dwa, int* ipdms,
	int* iinfo, float* rinfo,
	int* npcur, float* curcor, float* curvel,
	int* kxfo, float* rxfo, float* rhxfo,
	int* ixfo, int* iextf, int* icoord,
	int* nint, int* nrea, int* nsto, int* nstr,
	char chext[][MCHEXT], float state[][MSTATE],
	float gamma[][MGAMMA], float vellot[][MVELLOT],
	float* stor, int* ierr);

SOCKET connect_to_SAMS();