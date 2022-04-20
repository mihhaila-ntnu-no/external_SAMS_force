#include <math.h>
#include <stdbool.h>
#include <direct.h>
#include "simo_extfunc.h"
#include "csvwriter.h"
#include "getline.h"

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
	// Rename some input arguments for clarity
	int body_nr = iinfo[2]-1; // bodies are indexed starting with 1
	int nr_of_bodies = iinfo[4];
	int step_nr = iinfo[5];
	int nr_of_steps = iinfo[6];
	int substep_nr = iinfo[7];
	int nr_of_substeps = iinfo[8];
	int iteration_nr = iinfo[11];
	float time = rinfo[0];
	float dt = rinfo[1];
	float g = rinfo[2];
	
	// TODO test what happens if *ierr=-1. It's not guaranteed that the simulation will stop
	// TODO check that on this machine, float and double have a 2x size difference
	// The following arrays are declared as floats, but are in fact doubles spread over float pairs
	union
	{
		float single_precision[MGAMMA];
		double double_precision[MGAMMA / 2];
	} gamma_float;
	for (int i = 0; i < MGAMMA;i++)
		gamma_float.single_precision[i] = gamma[body_nr][i];
	// Because only half of gamma[] is available, the transformation matrix is irretrievable
	// Thus, the local velocities in vellot[] are in an unknown coordinate system
	union
	{
		float single_precision[MVELLOT];
		double double_precision[MVELLOT / 2];
	} vellot_float;
	for (int i = 0; i < MVELLOT;i++)
		vellot_float.single_precision[i] = vellot[body_nr][i];
	// First half of state[] contains displacements
	// The second half of state[] was supposed to contain the global velocities
	union
	{
		float single_precision[MSTATE];
		double double_precision[MSTATE / 2];
	} displacement_SIMA_t0; // timestep t=0
	for (int i = 0; i < MSTATE;i++)
		displacement_SIMA_t0.single_precision[i] = state[body_nr][i];

	// Declare memory structures for communication with SAMS
	union
	{
		double double_precision[1 + 6]; // time and central forces in 6DOF
		char character_array[(1 + 6) * sizeof(double)];
	} GFEXFO_to_SAMS;
	union
	{
		double double_precision[32];
		char character_array[32 * sizeof(double)];
	} SAMS_to_GFEXFO;

	// Static variables keep their value across timesteps and executions
	static SOCKET sams_tcp_socket;
	static int run_counter;
	static WSAPROTOCOL_INFO sams_tcp_socket_info;
	static CsvWriter* csv_writer;
	static double displacement_SIMA_tm1[6]; // Timestep t = -1
	static double displacement_SIMA_tm2[6]; // Timestep t = -2
	static double velocity_SAMS_tm1[6]; // Timestep t = -1
	static double mass_matrix_SIMA[6][6]; // [t=kg*10^3]
	static double mass_matrix_SAMS[6][6]; // [kg]
	static double stiffness_matrix_SIMA[6][6]; // Hydrostatic restoring force is modeled as spring stiffness
	static double stiffness_reference_SIMA[6]; // Hydrostatic equilibrium position

	// For now, assign constant forces regardless of the iteration or timestep
	if (step_nr < 10)
		stor[0] = 0.;
	else
		stor[0] = 10.; // [kN] surge force
	stor[1] = 0.;
	stor[2] = 0.;
	stor[3] = 0.;
	stor[4] = 0.;
	stor[5] = 0.;
	stor[6] = 0.;
	stor[7] = 0.;
	stor[8] = 0.;

	// Filenames and column names for logging
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
		"XGLB",
		"YGLB",
		"ZGLB",
		"FI",
		"THETA",
		"PSI",
		"velocity_SIMA_0",
		"velocity_SIMA_1",
		"velocity_SIMA_2",
		"velocity_SIMA_3",
		"velocity_SIMA_4",
		"velocity_SIMA_5",
		"acceleration_SIMA_0",
		"acceleration_SIMA_1",
		"acceleration_SIMA_2",
		"acceleration_SIMA_3",
		"acceleration_SIMA_4",
		"acceleration_SIMA_5",
		"inertial_force_SIMA_0",
		"inertial_force_SIMA_1",
		"inertial_force_SIMA_2",
		"inertial_force_SIMA_3",
		"inertial_force_SIMA_4",
		"inertial_force_SIMA_5",
		"hydrostatic_force_SIMA_0",
		"hydrostatic_force_SIMA_1",
		"hydrostatic_force_SIMA_2",
		"hydrostatic_force_SIMA_3",
		"hydrostatic_force_SIMA_4",
		"hydrostatic_force_SIMA_5",
		"SAMS_Time",
		"SAMS_positionX",
		"SAMS_positionY",
		"SAMS_positionZ",
		"SAMS_rotAxisX",
		"SAMS_rotAxisY",
		"SAMS_rotAxisZ",
		"SAMS_rotAngle",
		"SAMS_linVelX",
		"SAMS_linVelY",
		"SAMS_linVelZ",
		"SAMS_rotVelX",
		"SAMS_rotVelY",
		"SAMS_rotVelZ",
		"SAMS_linVelSurge",
		"SAMS_linVelSway",
		"SAMS_linVelHeave",
		"SAMS_rotVelRoll",
		"SAMS_rotVelPitch",
		"SAMS_rotVelYaw",
		"SAMS_totForceX",
		"SAMS_totForceY",
		"SAMS_totForceZ",
		"SAMS_totTorqueX",
		"SAMS_totTorqueY",
		"SAMS_totTorqueZ",
		"SAMS_totForceSurge",
		"SAMS_totForceSway",
		"SAMS_totForceHeave",
		"SAMS_totTorqueRoll",
		"SAMS_totTorquePitch",
		"SAMS_totTorqueYaw",
		"acceleration_SAMS_0",
		"acceleration_SAMS_1",
		"acceleration_SAMS_2",
		"acceleration_SAMS_3",
		"acceleration_SAMS_4",
		"acceleration_SAMS_5"
	};
	int nr_of_csv_ints = (int)sizeof(csv_ints_header) / sizeof(csv_ints_header[0]);
	int nr_of_csv_floats = (int)sizeof(csv_floats_header) / sizeof(csv_floats_header[0]);
	int nr_of_csv_doubles = (int)sizeof(csv_doubles_header) / sizeof(csv_doubles_header[0]);
	char* result_file_name = "gfexfo_sams.csv";
	
	char itconfig_file_path[MCHEXT];
	strncpy(itconfig_file_path, chext[0]+1, MCHEXT-1);
	char structure_file_path[MCHEXT];
	strncpy(structure_file_path, chext[1]+1, MCHEXT-1);
	
	run_counter++;

	// Missing kinematic state variables of current timestep
	double velocity_SIMA_t0[6];
	double velocity_SAMS_t0[6];
	double acceleration_SIMA[6];
	double acceleration_SAMS[6];
	double inertial_force_SIMA[6]; // [kN]
	double hydrostatic_force_SIMA[6]; // [kN]
	double radius_of_gyration_SAMS[3];
	double mass_matrix_difference[6][6];

	
	int iResult; // for error codes
	// First run, open the SAMS connection and write the log headers
	if (step_nr == 1 && substep_nr == 1 && iteration_nr == 0)
	{
		// TODO read the .itconfig, check that SAMS simulates one timestep less than SIMA
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

		char* line_buffer = NULL;
		size_t line_buffer_size = 0;
		int line_count = 0;
		ssize_t line_size;
		FILE* sys_dat_file = fopen("sys.dat", "r");
		if (!sys_dat_file)
		{
			printf("Error opening sys.dat\n");
			*ierr = 1;
			return;
		}
		do
		{
			getline(&line_buffer, &line_buffer_size, sys_dat_file);
		} while (strcmp(line_buffer, " MASS COEFFICIENTS\n"));
		getline(&line_buffer, &line_buffer_size, sys_dat_file); // Skip the separator line
		getline(&line_buffer, &line_buffer_size, sys_dat_file); // Skip the headers
		line_size = getline(&line_buffer, &line_buffer_size, sys_dat_file);
		if (line_size != 114)
			printf("Warning reading mass coefficients from sys.dat: unexpected length of line\n");
		iResult = sscanf
		(
			line_buffer,
			"%le%le%le%le%le%le%le\n",
			&mass_matrix_SIMA[0][0],
			&mass_matrix_SIMA[3][3],
			&mass_matrix_SIMA[4][3],
			&mass_matrix_SIMA[4][4],
			&mass_matrix_SIMA[5][3],
			&mass_matrix_SIMA[5][4],
			&mass_matrix_SIMA[5][5]
		);
		if (iResult != 7)
		{
			printf("Error parsing mass coefficients from sys.dat\n");
			*ierr = 1;
			return;
		}
		mass_matrix_SIMA[1][1] = mass_matrix_SIMA[0][0];
		mass_matrix_SIMA[2][2] = mass_matrix_SIMA[0][0];
		do
		{
			getline(&line_buffer, &line_buffer_size, sys_dat_file);
		} while (strcmp(line_buffer, " STIFFNESS REFERENCE\n"));
		getline(&line_buffer, &line_buffer_size, sys_dat_file); // Skip the separator line
		getline(&line_buffer, &line_buffer_size, sys_dat_file); // Skip the headers
		line_size = getline(&line_buffer, &line_buffer_size, sys_dat_file);
		if (line_size != 99)
			printf("Warning reading stiffness reference from sys.dat: unexpected length of line\n");
		iResult = sscanf
		(
			line_buffer,
			"%le%le%le%le%le%le\n",
			&stiffness_reference_SIMA[0],
			&stiffness_reference_SIMA[1],
			&stiffness_reference_SIMA[2],
			&stiffness_reference_SIMA[3],
			&stiffness_reference_SIMA[4],
			&stiffness_reference_SIMA[5]
		);
		if (iResult != 6)
		{
			printf("Error parsing stiffness coefficients from sys.dat\n");
			*ierr = 1;
			return;
		}
		do
		{
			getline(&line_buffer, &line_buffer_size, sys_dat_file);
		} while (strcmp(line_buffer, "'KMAT\n"));
		for (int i = 0;i < 6;i++)
		{
			line_size = getline(&line_buffer, &line_buffer_size, sys_dat_file);
			if (line_size != 86)
				printf("Warning reading stiffness matrix from sys.dat: unexpected length of line\n");
			iResult = sscanf
			(
				line_buffer,
				"%lf%lf%lf%lf%lf%lf\n",
				&stiffness_matrix_SIMA[i][0],
				&stiffness_matrix_SIMA[i][1],
				&stiffness_matrix_SIMA[i][2],
				&stiffness_matrix_SIMA[i][3],
				&stiffness_matrix_SIMA[i][4],
				&stiffness_matrix_SIMA[i][5]
			);
			if (iResult != 6)
			{
				printf("Error parsing stiffness coefficients from sys.dat\n");
				*ierr = 1;
				return;
			}
		}
		fclose(sys_dat_file);

		FILE* itconfig_file = fopen(itconfig_file_path, "r");
		if (!itconfig_file)
		{
			printf("Error opening:%s\n", itconfig_file_path);
			*ierr = 1;
			return;
		}
		fclose(itconfig_file);

		FILE* structure_file = fopen(structure_file_path, "r");
		if (!structure_file)
		{
			printf("Error opening:%s\n", structure_file_path);
			*ierr = 1;
			return;
		}
		do
		{
			getline(&line_buffer, &line_buffer_size, structure_file);
		} while (!strstr(line_buffer, "\"structureMass\""));
		iResult = sscanf(line_buffer, " \" structureMass \" : %lf ", &mass_matrix_SAMS[0][0]);
		if (iResult != 1)
		{
			// TODO show the actual file paths in these messages
			printf("Error parsing structure mass from structure file\n");
			*ierr = 1;
			return;
		}
		mass_matrix_SAMS[1][1] = mass_matrix_SAMS[0][0];
		mass_matrix_SAMS[2][2] = mass_matrix_SAMS[0][0];

		// Do not know whether the mass or the radii come first, so search from the top
		fclose(structure_file);
		structure_file = fopen(structure_file_path, "r");
		do
		{
			getline(&line_buffer, &line_buffer_size, structure_file);
		} while (!strstr(line_buffer, "\"structureRadiusOfGyration\""));
		iResult = sscanf
		(
			line_buffer,
			" \" structureRadiusOfGyration \" : [ %lf , %lf , %lf",
			&radius_of_gyration_SAMS[0],
			&radius_of_gyration_SAMS[1],
			&radius_of_gyration_SAMS[2]
		);
		if (iResult != 3)
		{
			// TODO show the actual file paths in these messages
			printf("Error parsing gyration radii from structure file\n");
			*ierr = 1;
			return;
		}
		fclose(structure_file);
		
		free(line_buffer);
		// Calculate the mass moments of inertia from the radii of gyration
		for (int i = 0;i < 3;i++)
			mass_matrix_SAMS[3 + i][3 + i] = pow(radius_of_gyration_SAMS[i], 2) * mass_matrix_SAMS[i][i];
		for (int i = 0;i < 6;i++)
		{
			for (int j = 0;j < 6;j++)
			{
				if (mass_matrix_SAMS[i][j] && mass_matrix_SIMA[i][j] )
				{
					mass_matrix_difference[i][j] =
						(mass_matrix_SIMA[i][j] * 1000 - mass_matrix_SAMS[i][j])
						/ mass_matrix_SIMA[i][j];
					if (mass_matrix_difference[i][j] > 0.01)
					{
						printf
						(
							"At position [%d][%d], the mass matrix has %lf in SIMA and %lf in SAMS, relative difference %le\n",
							i,
							j,
							mass_matrix_SIMA[i][j] * 1000,
							mass_matrix_SAMS[i][j],
							mass_matrix_difference[i][j]
						);
						*ierr = 1;
						return;
					}
				}
			}
		}
	}
	
	// Calculate what we can in 6DoF without having received anything from SAMS
	for (int i = 0;i < 6;i++)
	{
		// Back-initialize the SIMA displacements to avoid initial jerk
		if (run_counter == 1)
		{
			displacement_SIMA_tm1[i] = displacement_SIMA_t0.double_precision[i];
			displacement_SIMA_tm2[i] = displacement_SIMA_t0.double_precision[i];
		}
		// Calculate velocities and accelerations as backward finite differences of displacement.
		// In the simulation beginning, unknown values are assumed 0
		velocity_SIMA_t0[i] =
			(
				-displacement_SIMA_tm1[i]
				+ displacement_SIMA_t0.double_precision[i]
				) / dt;
		acceleration_SIMA[i] =
			(
				displacement_SIMA_tm2[i]
				- 2 * displacement_SIMA_tm1[i]
				+ displacement_SIMA_t0.double_precision[i]
				) / pow(dt, 2);
		// Calculate inertial forces
		inertial_force_SIMA[i] = acceleration_SIMA[i] * mass_matrix_SIMA[i][i];
		hydrostatic_force_SIMA[i] = (stiffness_reference_SIMA[i] - displacement_SIMA_t0.double_precision[i]) * stiffness_matrix_SIMA[i][i];
	}
	// TODO Test resilience for extra timesteps and iterations. Right now, both are disabled through SIMA settings
	if (iteration_nr == 0)
	{
		if (sams_tcp_socket == INVALID_SOCKET)
		{
			printf("Invalid socket. Returning default value 1\n");
			*ierr = 1;
			return;
		}

		// Receive kinematic and dynamic information from SAMS
		iResult = recv(sams_tcp_socket, SAMS_to_GFEXFO.character_array, sizeof(SAMS_to_GFEXFO.double_precision), 0);
		if (iResult == 0)
		{
			// This never seems to happen
			printf("Connection closed.\n");
			*ierr = 1;
			return;
		}
		else if (iResult < 0)
		{
			printf("Receive failed: %d\n", WSAGetLastError());
			*ierr = 1;
			return;
		}

		for (int i = 0;i < 6;i++)
		{
			// Rename some received variables for clarity
			velocity_SAMS_t0[i] = SAMS_to_GFEXFO.double_precision[8 + i];
			// Compute  accelerations as backward finite differences of velocity.
			acceleration_SAMS[i] =
				(
					-velocity_SAMS_tm1[i]
					+ velocity_SAMS_t0[i]
					) / dt;
		}

		// SIMO receives forces through stor[] in kN, SAMS through TCP in N
		GFEXFO_to_SAMS.double_precision[0] = (double)time; // [s]
		GFEXFO_to_SAMS.double_precision[1] = stor[0] * 1000;
		GFEXFO_to_SAMS.double_precision[2] = stor[1] * 1000;
		GFEXFO_to_SAMS.double_precision[3] = stor[2] * 1000;
		GFEXFO_to_SAMS.double_precision[4] = stor[3] * 1000;
		GFEXFO_to_SAMS.double_precision[5] = stor[4] * 1000;
		GFEXFO_to_SAMS.double_precision[6] = stor[5] * 1000;
		iResult = send(sams_tcp_socket, GFEXFO_to_SAMS.character_array, sizeof(GFEXFO_to_SAMS.double_precision), 0);
		if (iResult == SOCKET_ERROR)
		{
			printf("Send failed: %d\n", WSAGetLastError());
			closesocket(sams_tcp_socket);
			WSACleanup();
			*ierr = 1;
			return;
		}

		// Log the variables for this timestep
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
			displacement_SIMA_t0.double_precision[0],
			displacement_SIMA_t0.double_precision[1],
			displacement_SIMA_t0.double_precision[2],
			displacement_SIMA_t0.double_precision[3],
			displacement_SIMA_t0.double_precision[4],
			displacement_SIMA_t0.double_precision[5],
			velocity_SIMA_t0[0],
			velocity_SIMA_t0[1],
			velocity_SIMA_t0[2],
			velocity_SIMA_t0[3],
			velocity_SIMA_t0[4],
			velocity_SIMA_t0[5],
			acceleration_SIMA[0],
			acceleration_SIMA[1],
			acceleration_SIMA[2],
			acceleration_SIMA[3],
			acceleration_SIMA[4],
			acceleration_SIMA[5],
			inertial_force_SIMA[0],
			inertial_force_SIMA[1],
			inertial_force_SIMA[2],
			inertial_force_SIMA[3],
			inertial_force_SIMA[4],
			inertial_force_SIMA[5],
			hydrostatic_force_SIMA[0],
			hydrostatic_force_SIMA[1],
			hydrostatic_force_SIMA[2],
			hydrostatic_force_SIMA[3],
			hydrostatic_force_SIMA[4],
			hydrostatic_force_SIMA[5],
			SAMS_to_GFEXFO.double_precision[0], // SAMS_Time
			SAMS_to_GFEXFO.double_precision[1], // SAMS_positionX
			SAMS_to_GFEXFO.double_precision[2], // SAMS_positionY
			SAMS_to_GFEXFO.double_precision[3], // SAMS_positionZ
			SAMS_to_GFEXFO.double_precision[4], // SAMS_rotAxisX
			SAMS_to_GFEXFO.double_precision[5], // SAMS_rotAxisY
			SAMS_to_GFEXFO.double_precision[6], // SAMS_rotAxisZ
			SAMS_to_GFEXFO.double_precision[7], // SAMS_rotAngle
			SAMS_to_GFEXFO.double_precision[8], // SAMS_linVelX
			SAMS_to_GFEXFO.double_precision[9], // SAMS_linVelY
			SAMS_to_GFEXFO.double_precision[10], // SAMS_linVelZ
			SAMS_to_GFEXFO.double_precision[11], // SAMS_rotVelX
			SAMS_to_GFEXFO.double_precision[12], // SAMS_rotVelY
			SAMS_to_GFEXFO.double_precision[13], // SAMS_rotVelZ
			SAMS_to_GFEXFO.double_precision[14], // SAMS_linVelSurge
			SAMS_to_GFEXFO.double_precision[15], // SAMS_linVelSway
			SAMS_to_GFEXFO.double_precision[16], // SAMS_linVelHeave
			SAMS_to_GFEXFO.double_precision[17], // SAMS_rotVelRoll
			SAMS_to_GFEXFO.double_precision[18], // SAMS_rotVelPitch
			SAMS_to_GFEXFO.double_precision[19], // SAMS_rotVelYaw
			SAMS_to_GFEXFO.double_precision[20], // SAMS_totForceX
			SAMS_to_GFEXFO.double_precision[21], // SAMS_totForceY
			SAMS_to_GFEXFO.double_precision[22], // SAMS_totForceZ
			SAMS_to_GFEXFO.double_precision[23], // SAMS_totTorqueX
			SAMS_to_GFEXFO.double_precision[24], // SAMS_totTorqueY
			SAMS_to_GFEXFO.double_precision[25], // SAMS_totTorqueZ
			SAMS_to_GFEXFO.double_precision[26], // SAMS_totForceSurge
			SAMS_to_GFEXFO.double_precision[27], // SAMS_totForceSway
			SAMS_to_GFEXFO.double_precision[28], // SAMS_totForceHeave
			SAMS_to_GFEXFO.double_precision[29], // SAMS_totForceRoll
			SAMS_to_GFEXFO.double_precision[30], // SAMS_totForcePitch
			SAMS_to_GFEXFO.double_precision[31], // SAMS_totForceYaw
			acceleration_SAMS[0],
			acceleration_SAMS[1],
			acceleration_SAMS[2],
			acceleration_SAMS[3],
			acceleration_SAMS[4],
			acceleration_SAMS[5]
		};
		CsvWriter_nextRow(csv_writer);
		char log_buffer[50];
		for (int i = 0; i < nr_of_csv_ints; i++)
		{
			sprintf(log_buffer, "%d", csv_ints[i]);
			if (CsvWriter_writeField(csv_writer, log_buffer))
			{
				printf("Error writing .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
				*ierr = 1;
				return;
			}
		}
		for (int i = 0; i < nr_of_csv_floats; i++)
		{
			sprintf(log_buffer, "%f", csv_floats[i]);
			if (CsvWriter_writeField(csv_writer, log_buffer))
			{
				printf("Error writing .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
				*ierr = 1;
				return;
			}
		}
		for (int i = 0; i < nr_of_csv_doubles; i++)
		{
			sprintf(log_buffer, "%f", csv_doubles[i]);
			if (CsvWriter_writeField(csv_writer, log_buffer))
			{
				printf("Error writing .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
				*ierr = 1;
				return;
			}
		}

		// Calculations done, prepare for the next timestep
		for (int i = 0; i < 6; i++)
		{
			displacement_SIMA_tm2[i] = displacement_SIMA_tm1[i];
			displacement_SIMA_tm1[i] = displacement_SIMA_t0.double_precision[i];
			velocity_SAMS_tm1[i] = velocity_SAMS_t0[i];
		}
	}
	// Last run, close the connection and logs
	if (step_nr == nr_of_steps && substep_nr == nr_of_substeps && iteration_nr == 0)
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
