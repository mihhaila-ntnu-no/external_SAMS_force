#include <math.h>
#include <stdbool.h>
#include <direct.h>
#include "simo_extfunc.h"
#include "csvwriter.h"

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
	int body_nr = iinfo[2]-1; // bodies are indexed starting with 1
	int nr_of_bodies = iinfo[4];
	int step_nr = iinfo[5];
	int substep_nr = iinfo[7];
	int iteration_nr = iinfo[11];

	static SOCKET sams_tcp_socket;
	static int run_counter;
	static WSAPROTOCOL_INFO sams_tcp_socket_info;
	static CsvWriter* csv_writer;
	static FILE* log_file;

	// For now, assign constant forces regardless of the iteration or timestep
	stor[0] = 10.; // [kN] surge force
	stor[1] = 0.;
	stor[2] = 0.;
	stor[3] = 0.;
	stor[4] = 0.;
	stor[5] = 0.;
	stor[6] = 0.;
	stor[7] = 0.;
	stor[8] = 0.;

	// TODO check that on this machine, float and double have a 2x size difference
	// Some single-precision float pairs in the argument are actually double-precision floats!
	union
	{
		float single_precision[MSTATE];
		double double_precision[MSTATE / 2];
	} state_float;
	for (int i = 0; i < MSTATE;i++)
	{
		state_float.single_precision[i] = state[body_nr][i];
	}
	union
	{
		float single_precision[MVELLOT];
		double double_precision[MVELLOT / 2];
	} vellot_float;
	for (int i = 0; i < MVELLOT;i++)
	{
		vellot_float.single_precision[i] = vellot[body_nr][i];
	}
	union
	{
		float single_precision[MGAMMA];
		double double_precision[MGAMMA / 2];
	} gamma_float;
	for (int i = 0; i < MGAMMA;i++)
	{
		gamma_float.single_precision[i] = gamma[body_nr][i];
	}

	char* csv_ints_header[] =
	{
		"iwa",
		"ipdms",
		"IMODE",
		"MODULE",
		"IBDY",
		"IBDTYP",
		"NBDY",
		"ISTEP",
		"NSTEP",
		"IEXTRA",
		"NEXTRA",
		"IGRAV",
		"ISTORE",
		"ITER",
		"npcur",
		"kxfo",
		"ixfo",
		"iextf",
		"icoord",
		"nint",
		"nrea",
		"nsto",
		"nstr",
	};
	char* csv_floats_header[] =
	{
		"rwa",
		"TIME",
		"DT",
		"GRAV",
		"RHOW",
		"RHOA",
		"DEPTH",
		"curcor",
		"curvel",
		"rxfo",
		"rhxfo",
		"stor_0",
		"stor_1",
		"stor_2",
		"stor_3",
		"stor_4",
		"stor_5",
		"stor_6",
		"stor_7",
		"stor_8"
	};
	char* csv_doubles_header[] =
	{
		"dwa",
		"state_0",
		"state_1",
		"state_2",
		"state_3",
		"state_4",
		"state_5",
		"gamma_0",
		"gamma_1",
		"gamma_2",
		"gamma_3",
		"vellot_0",
		"vellot_1",
		"vellot_2",
		"SAMS_time",
		"SAMS_force_X",
		"SAMS_force_Y",
		"SAMS_force_Z",
		"SAMS_moment_X",
		"SAMS_moment_Y",
		"SAMS_moment_Z"
	};
	char* result_file_name = "gfexfo_sams.csv";
	char* log_file_name = "gfexfo_sams.log";


	// Connect numeric parameters into typed arrays
	int csv_ints[] =
	{
		*iwa,
		*ipdms,
		iinfo[0],
		iinfo[1],
		iinfo[2],
		iinfo[3],
		iinfo[4],
		iinfo[5],
		iinfo[6],
		iinfo[7],
		iinfo[8],
		iinfo[9],
		iinfo[10],
		iinfo[11],
		*npcur,
		*kxfo,
		*ixfo,
		*iextf,
		*icoord,
		*nint,
		*nrea,
		*nsto,
		*nstr,
	};
	float csv_floats[] =
	{
		*rwa,
		rinfo[0],
		rinfo[1],
		rinfo[2],
		rinfo[3],
		rinfo[4],
		rinfo[5],
		*curcor,
		*curvel,
		*rxfo,
		*rhxfo,
		stor[0],
		stor[1],
		stor[2],
		stor[3],
		stor[4],
		stor[5],
		stor[6],
		stor[7],
		stor[8]
	};
	double csv_doubles[] =
	{
		*dwa,
		state_float.double_precision[0],
		state_float.double_precision[1],
		state_float.double_precision[2],
		state_float.double_precision[3],
		state_float.double_precision[4],
		state_float.double_precision[5],
		gamma_float.double_precision[0],
		gamma_float.double_precision[1],
		gamma_float.double_precision[2],
		gamma_float.double_precision[3],
		vellot_float.double_precision[0],
		vellot_float.double_precision[1],
		vellot_float.double_precision[2]
	};
	int nr_of_csv_ints = (int)sizeof(csv_ints) / sizeof(csv_ints[0]);
	int nr_of_csv_floats = (int)sizeof(csv_floats) / sizeof(csv_floats[0]);
	int nr_of_csv_doubles = (int)sizeof(csv_doubles) / sizeof(csv_doubles[0]);


	run_counter++;

	int iResult; // for error codes

	if (step_nr == 1 && substep_nr == 1 && iteration_nr == 0)
	{
		log_file = fopen(log_file_name, "a+");
		sams_tcp_socket = connect_to_SAMS();
		csv_writer = CsvWriter_new(result_file_name, ";", 0);
		for (int i = 0; i < nr_of_csv_ints; i++)
		{
			if (CsvWriter_writeField(csv_writer, csv_ints_header[i]))
			{
				printf("Error writing .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
				*ierr = 1;
				return;
			}
		}
		for (int i = 0; i < nr_of_csv_floats; i++)
		{
			if (CsvWriter_writeField(csv_writer, csv_floats_header[i]))
			{
				printf("Error writing .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
				*ierr = 1;
				return;
			}
		}
		for (int i = 0; i < nr_of_csv_doubles; i++)
		{
			if (CsvWriter_writeField(csv_writer, csv_doubles_header[i]))
			{
				printf("Error writing .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
				*ierr = 1;
				return;
			}
		}
	}

	if (iteration_nr == 0)
	{
		CsvWriter_nextRow(csv_writer);

		// Collect numeric parameters into string arrays
		// TODO don't need to collect an array of string buffers, could just reuse one every loop
		char csv_ints_printed[23][17];
		char csv_floats_printed[20][50]; // 16 recommended + 1 terminating
		char csv_doubles_printed[14][50]; // 24 recommended
		for (int i = 0; i < nr_of_csv_ints; i++)
		{
			sprintf(csv_ints_printed[i], "%d", csv_ints[i]);
		}
		for (int i = 0; i < nr_of_csv_floats; i++)
		{
			sprintf(csv_floats_printed[i], "%f", csv_floats[i]);
		}
		for (int i = 0; i < nr_of_csv_doubles; i++)
		{
			sprintf(csv_doubles_printed[i], "%f", csv_doubles[i]);
		}
		// Print row of values
		for (int i = 0; i < nr_of_csv_ints; i++)
		{
			if (CsvWriter_writeField(csv_writer, csv_ints_printed[i]))
			{
				printf("Error writing .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
				*ierr = 1;
				return;
			}
		}
		for (int i = 0; i < nr_of_csv_floats; i++)
		{
			if (CsvWriter_writeField(csv_writer, csv_floats_printed[i]))
			{
				printf("Error writing .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
				*ierr = 1;
				return;
			}
		}
		for (int i = 0; i < nr_of_csv_doubles; i++)
		{
			if (CsvWriter_writeField(csv_writer, csv_doubles_printed[i]))
			{
				printf("Error writing .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
				*ierr = 1;
				return;
			}
		}

		if (sams_tcp_socket == INVALID_SOCKET)
		{
			printf("Invalid socket. Returning default value 1\n");
			*ierr = 1;
			return;
		}
		union
		{
			double double_precision[1+6]; // time and central forces in 6DOF
			char character_array[(1+6) * sizeof(double)];
		} input_to_SAMS;
		input_to_SAMS.double_precision[0] = (double)rinfo[0]; // [s] time
		input_to_SAMS.double_precision[1] = 10000.; // [N] surge force
		input_to_SAMS.double_precision[2] = 0.;
		input_to_SAMS.double_precision[3] = 0.;
		input_to_SAMS.double_precision[4] = 0.;
		input_to_SAMS.double_precision[5] = 0.;
		input_to_SAMS.double_precision[6] = 0.;
		iResult = send(sams_tcp_socket, input_to_SAMS.character_array, sizeof(input_to_SAMS.double_precision), 0);
		if (iResult == SOCKET_ERROR)
		{
			printf("Send failed: %d\n", WSAGetLastError());
			closesocket(sams_tcp_socket);
			WSACleanup();
			*ierr = 1;
			return;
		}

		union
		{
			double double_precision[32];
			char character_array[32 * sizeof(double)];
		} output_from_SAMS;
		iResult = recv(sams_tcp_socket, output_from_SAMS.character_array, sizeof(output_from_SAMS.double_precision), 0);

		if (iResult > 0)
		{
			// use output_from_SAMS.double_precision
		}
		else if (iResult == 0)
		{
			printf("Connection closed.\n");
		}
		else
		{
			printf("Receive failed: %d\n", WSAGetLastError());
			*ierr = 1;
			return;
		}
	}
	if (step_nr == iinfo[6] && substep_nr == iinfo[8] && iteration_nr == 0)
	{
		iResult = WSACleanup();
		if (iResult == SOCKET_ERROR)
		{
			printf("Error closing connection: %d\n", WSAGetLastError());
			*ierr = 1;
			return;
		}
		CsvWriter_destroy(csv_writer);
		printf("Disconnected from SAMS, results saved in %s\n",result_file_name);
		fclose(log_file); // TODO close log file even in case of error and return
	}


	*ierr = 0;
}

SOCKET connect_to_SAMS()
{
	WSADATA wsaData;
	int iResult = WSAStartup(MAKEWORD(2, 2), &wsaData);
	if (iResult != 0)
	{
		printf("WSAStartup failed: %d\n", iResult);
		return INVALID_SOCKET;
	}
	struct addrinfo* result = NULL,
		* ptr = NULL,
		hints;
	ZeroMemory(&hints, sizeof(hints));
	hints.ai_family = AF_INET;
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_protocol = IPPROTO_TCP;
	iResult = getaddrinfo("localhost", "27015", &hints, &result);
	if (iResult != 0)
	{
		printf("getaddrinfo failed: %d\n", iResult);
		WSACleanup();
		return INVALID_SOCKET;
	}
	SOCKET ConnectSocket = INVALID_SOCKET;
	ptr = result;
	bool connection_found = 0;
	while (ptr != NULL)
	{
		ConnectSocket = socket(ptr->ai_family, ptr->ai_socktype, ptr->ai_protocol);
		if (ConnectSocket == INVALID_SOCKET)
		{
			printf("Error at socket(): %ld\n", WSAGetLastError());
			ptr = ptr->ai_next;
			continue;
		}
		iResult = connect(ConnectSocket, ptr->ai_addr, (int)ptr->ai_addrlen);
		if (iResult == SOCKET_ERROR)
		{
			printf("Socket error: %ld\n", WSAGetLastError());
			ptr = ptr->ai_next;
			continue;
		}
		connection_found = 1;
		break;
	}
	freeaddrinfo(result);
	if (connection_found)
	{
		printf("Connected to SAMS\n");
		return ConnectSocket;
	}
	else
	{
		printf("Unable to connect to server\n");
		WSACleanup();
		return INVALID_SOCKET;
	}
}