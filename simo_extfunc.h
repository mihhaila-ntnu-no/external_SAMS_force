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

int send_to_SAMS(SOCKET sams_tcp_socket, double time, double central_forces[6]);

int receive_from_SAMS(SOCKET sams_tcp_socket, double* SAMS_time, double displacement_SAMS[6], double SIMA_global_to_SAMS_body[3][3]);

int read_from_SAMS(char* SAMS_resultfile_path, int substep_nr_overall, int nr_of_substeps_overall, double* SAMS_txt_time, double displacement_SAMS_t0[6], double velocity_SAMS_t0[6], double F_ice_global[6]);

int get_matrix_from_axis_angle(double e1, double e2, double e3, double theta, double SAMS_global_to_SAMS_body[3][3]);

int get_joint_rotations_from_matrix(double SAMS_global_to_SAMS_body[3][3], double joint_rotation[3]);