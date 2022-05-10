#define _USE_MATH_DEFINES
#include <math.h>
#include <stdbool.h>
#include <direct.h>
#include "simo_extfunc.h"
#include "csvwriter.h"
#include "getline.h"

void CAL_CONV gfexfo_
(
	int* iwa,
	float* rwa,
	double* dwa,
	int* ipdms,
	int* iinfo,
	float* rinfo,
	int* npcur,
	float* curcor,
	float* curvel,
	int* kxfo,
	float* rxfo,
	float* rhxfo,
	int* ixfo,
	int* iextf,
	int* icoord,
	int* nint,
	int* nrea,
	int* nsto,
	int* nstr,
	char chext[][MCHEXT],
	float state[][MSTATE],
	float gamma[][MGAMMA],
	float vellot[][MVELLOT],
	float* stor,
	int* ierr
)
{
	// *ierr will be incremented in case of a warning, so at the end of execution the warnings can be counted
	// In case of an error, *ierr will be set to a negative value and gfexfo_() will return
	// The error values decrease sequentially as they appear in the code, but this hasn't shown up in SIMA so far
	*ierr = 0;
	// If this is not the first iteration of the timestep, return the forces from the first one
	if (iinfo[11] > 0) // iinfo[11] = iteration nr, starting from 0
	{
		// TODO return SAMS_ice_force from a previous run
		stor[0] = 0.; // [kN] surge
		stor[1] = 0.; // [kN] sway
		stor[2] = 0.; // [kN] heave
		stor[3] = 0.; // [MN*m] roll
		stor[4] = 0.; // [MN*m] pitch
		stor[5] = 0.; // [MN*m] yaw
		stor[6] = 0.; // reserved, don't touch
		stor[7] = 0.; // reserved, don't touch
		stor[8] = 0.; // reserved, don't touch
		return;
	}
	// Function argument aliases
	int body_nr = iinfo[2]-1; // bodies are indexed starting with 1
	int step_nr = iinfo[5];
	int nr_of_steps = iinfo[6];
	int substep_nr = iinfo[7];
	int nr_of_substeps = iinfo[8];
	float time = rinfo[0];
	float dt = rinfo[1];
	
	// TODO check that on this machine, float and double have a 2x size difference
	// TODO determine what should be cleaned up, do it on every error. WSA, buffer free, etc
	
	// The following arrays are declared as (single-precision) floats in the function definition
	// In fact, they contain doubles spread over float pairs
	union
	{
		float single_precision[MGAMMA];
		double double_precision[MGAMMA / 2]; // 9/2=4 for integers
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
	} displacement_SIMA_t0; // timestep t=0. Rotations in radians
	for (int i = 0; i < MSTATE;i++)
		displacement_SIMA_t0.single_precision[i] = state[body_nr][i];

	// Static variables keep their value across timesteps and executions
	static int run_counter;
	run_counter++;
	static SOCKET sams_tcp_socket;
	static WSAPROTOCOL_INFO sams_tcp_socket_info;
	static CsvWriter* csv_writer;
	static double displacement_SIMA_tm1[6]; // Timestep t = -1.
	static double displacement_SIMA_tm2[6]; // Timestep t = -2
	static double velocity_SAMS_tm1[6]; // Timestep t = -1
	static double mass_matrix_SIMA[6][6]; // [t=kg*10^3]
	static double mass_matrix_SAMS[6][6]; // [kg]
	static double stiffness_matrix_SIMA[6][6]; // Hydrostatic restoring force is modeled as spring stiffness
	static double stiffness_reference_SIMA[6]; // Hydrostatic equilibrium position
	static char SAMS_resultfile_path[MCHEXT];
	static FILE* SAMS_result_file;
	static double SAMS_ice_force[6];

	// Assign forces external to SIMA
	{
		stor[0] = 0.; // [kN] surge
		stor[1] = 0.; // [kN] sway
		stor[2] = 0.; // [kN] heave
		stor[3] = 0.; // [MN*m] roll
		stor[4] = 0.; // [MN*m] pitch
		stor[5] = 0.; // [MN*m] yaw
		stor[6] = 0.; // reserved, don't touch
		stor[7] = 0.; // reserved, don't touch
		stor[8] = 0.; // reserved, don't touch
	}

	// Log file and column names
	char* gfexfo_result_file_name = "gfexfo_sams.csv";
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
		"displacement_SAMS_0",
		"displacement_SAMS_1",
		"displacement_SAMS_2",
		"displacement_SAMS_3",
		"displacement_SAMS_4",
		"displacement_SAMS_5",
		"acceleration_SAMS_0",
		"acceleration_SAMS_1",
		"acceleration_SAMS_2",
		"acceleration_SAMS_3",
		"acceleration_SAMS_4",
		"acceleration_SAMS_5",
		"inertial_force_SAMS_0",
		"inertial_force_SAMS_1",
		"inertial_force_SAMS_2",
		"inertial_force_SAMS_3",
		"inertial_force_SAMS_4",
		"inertial_force_SAMS_5",
		"hydrostatic_force_SAMS_0",
		"hydrostatic_force_SAMS_1",
		"hydrostatic_force_SAMS_2",
		"hydrostatic_force_SAMS_3",
		"hydrostatic_force_SAMS_4",
		"hydrostatic_force_SAMS_5",
		"rotation_matrix_SIMA_0_0",
		"rotation_matrix_SIMA_0_1",
		"rotation_matrix_SIMA_0_2",
		"rotation_matrix_SIMA_1_0",
		"rotation_matrix_SIMA_1_1",
		"rotation_matrix_SIMA_1_2",
		"rotation_matrix_SIMA_2_0",
		"rotation_matrix_SIMA_2_1",
		"rotation_matrix_SIMA_2_2",
		"gamma_0",
		"gamma_1",
		"gamma_2",
		"gamma_3",
		"rotation_matrix_SAMS_0_0",
		"rotation_matrix_SAMS_0_1",
		"rotation_matrix_SAMS_0_2",
		"rotation_matrix_SAMS_1_0",
		"rotation_matrix_SAMS_1_1",
		"rotation_matrix_SAMS_1_2",
		"rotation_matrix_SAMS_2_0",
		"rotation_matrix_SAMS_2_1",
		"rotation_matrix_SAMS_2_2",
		"SAMS_ice_force_0",
		"SAMS_ice_force_1",
		"SAMS_ice_force_2",
		"SAMS_ice_force_3",
		"SAMS_ice_force_4",
		"SAMS_ice_force_5",
	};
	int nr_of_csv_ints = (int)sizeof(csv_ints_header) / sizeof(csv_ints_header[0]);
	int nr_of_csv_floats = (int)sizeof(csv_floats_header) / sizeof(csv_floats_header[0]);
	int nr_of_csv_doubles = (int)sizeof(csv_doubles_header) / sizeof(csv_doubles_header[0]);

	// Working variables
	double velocity_SIMA_t0[6];
	double velocity_SAMS_t0[6];
	double acceleration_SIMA[6];
	double displacement_SAMS[6];
	double acceleration_SAMS[6];
	double inertial_force_SIMA[6]; // [kN, MN*m]
	double inertial_force_SAMS[6]; // [N]
	double hydrostatic_force_SIMA[6]; // [kN, MN*m]
	double hydrostatic_force_SAMS[6]; // [N]
	double rotation_matrix_SIMA[3][3];
	double rotation_matrix_SAMS[3][3];
	double SAMS_txt_output[47];
	double F_coupled_global[6]; // Sea forces on structure, as opposed to ice forces
	double F_coupled_local[6];
	double SAMS_time;
	int iResult; // For error codes
	char* line_buffer = NULL; // for getline()
	size_t line_buffer_size = 0; // for getline()
	
	// Calculate kinematic state, state-dependent forces, global->body transformation matrix
	{
		// Back-initialize the SIMA displacements to avoid calculating an unphysical initial jerk
		if (run_counter == 1)
		{
			for (int i = 0;i < 6;i++)
			{
				displacement_SIMA_tm1[i] = displacement_SIMA_t0.double_precision[i];
				displacement_SIMA_tm2[i] = displacement_SIMA_t0.double_precision[i];
			}
		}
		// Calculate kinematic state in global coordinates from displacement history
		for (int i = 0;i < 3;i++)
		{
			double rotation_velocity_negative_jump;
			double rotation_velocity_positive_jump;
			double rotation_acceleration_negative_jump;
			double rotation_acceleration_positive_jump;
			// Calculate translational velocities and accelerations as backward finite differences of displacement.
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
			// Calculate rotational velocities and accelerations, check for 0-360 = 0-2pi transition jump
			velocity_SIMA_t0[3 + i] =
				(
					-displacement_SIMA_tm1[3 + i]
					+ displacement_SIMA_t0.double_precision[3 + i]
					) / dt;
			rotation_velocity_negative_jump =
				(
					-displacement_SIMA_tm1[3 + i]
					+ displacement_SIMA_t0.double_precision[3 + i]
					- 2 * M_PI
					) / dt;
			rotation_velocity_positive_jump =
				(
					-displacement_SIMA_tm1[3 + i]
					+ displacement_SIMA_t0.double_precision[3 + i]
					+ 2 * M_PI
					) / dt;
			if (abs(velocity_SIMA_t0[3 + i]) > abs(rotation_velocity_negative_jump))
				velocity_SIMA_t0[3 + i] = rotation_velocity_negative_jump;
			if (abs(velocity_SIMA_t0[3 + i]) > abs(rotation_velocity_positive_jump))
				velocity_SIMA_t0[3 + i] = rotation_velocity_positive_jump;
			acceleration_SIMA[3 + i] =
				(
					displacement_SIMA_tm2[3 + i]
					- 2 * displacement_SIMA_tm1[3 + i]
					+ displacement_SIMA_t0.double_precision[3 + i]
					) / pow(dt, 2);
			rotation_acceleration_negative_jump =
				(
					displacement_SIMA_tm2[3 + i]
					- 2 * displacement_SIMA_tm1[3 + i]
					+ displacement_SIMA_t0.double_precision[3 + i]
					- 2 * M_PI
					) / pow(dt, 2);
			rotation_acceleration_positive_jump =
				(
					displacement_SIMA_tm2[3 + i]
					- 2 * displacement_SIMA_tm1[3 + i]
					+ displacement_SIMA_t0.double_precision[3 + i]
					+ 2 * M_PI
					) / pow(dt, 2);
			if (abs(acceleration_SIMA[3 + i]) > abs(rotation_acceleration_negative_jump))
				acceleration_SIMA[3 + i] = rotation_acceleration_negative_jump;
			if (abs(acceleration_SIMA[3 + i]) > abs(rotation_acceleration_positive_jump))
				acceleration_SIMA[3 + i] = rotation_acceleration_positive_jump;

		}
		// Rename the sequential angular displacements for clarity
		double phi = displacement_SIMA_t0.double_precision[3];
		double theta = displacement_SIMA_t0.double_precision[4];
		double psi = displacement_SIMA_t0.double_precision[5];
		// Calculate the global -> body coordinate rotation matrix
		rotation_matrix_SIMA[0][0] = cos(psi) * cos(theta);
		rotation_matrix_SIMA[0][1] = -sin(psi) * cos(phi) + cos(psi) * sin(theta) * sin(phi);
		rotation_matrix_SIMA[0][2] = sin(psi) * sin(phi) + cos(psi) * sin(theta) * cos(phi);
		rotation_matrix_SIMA[1][0] = sin(psi) * cos(theta);
		rotation_matrix_SIMA[1][1] = cos(psi) * cos(phi) + sin(psi) * sin(theta) * sin(phi);
		rotation_matrix_SIMA[1][2] = -cos(psi) * sin(phi) + sin(psi) * sin(theta) * cos(phi);
		rotation_matrix_SIMA[2][0] = -sin(theta);
		rotation_matrix_SIMA[2][1] = cos(theta) * sin(phi);
		rotation_matrix_SIMA[2][2] = cos(theta) * cos(phi);
		// Calculate forces acting on structure, in global coordinates
		for (int i = 0;i < 6;i++)
		{
			inertial_force_SIMA[i] = acceleration_SIMA[i] * mass_matrix_SIMA[i][i];
			hydrostatic_force_SIMA[i] =
				(stiffness_reference_SIMA[i] - displacement_SIMA_t0.double_precision[i])
				* stiffness_matrix_SIMA[i][i];
			F_coupled_global[i] = inertial_force_SIMA[i] + hydrostatic_force_SIMA[i];
		}
		// Transform the forces from global coordinate system to local
		for (int i = 0;i < 3;i++)
		{
			F_coupled_local[i] = 0;
			for (int j = 0;j < 3;j++)
				F_coupled_local[i] += rotation_matrix_SIMA[j][i] * F_coupled_global[j];
			// The rotational forces are the same in both coordinate systems
			F_coupled_local[3 + i] = F_coupled_global[3 + i];
		}
	}

	// First run
	if (run_counter == 1)
	{
		sams_tcp_socket = connect_to_SAMS();
		if (sams_tcp_socket == INVALID_SOCKET)
		{
			printf("Unable to connect to SAMS\n");
			WSACleanup();
			*ierr = -27;
			return;
		}
		// The initial timestep i=0 happened without calling gfexfo_(), nothing was logged from that timestep
		// Similarly, the initial timestep i=0 is received from SAMS, but nothing is done with it
		// Instead, it will be shortly overwritten by the information of i=1
		iResult = receive_from_SAMS(sams_tcp_socket, &SAMS_time, displacement_SAMS, velocity_SAMS_t0, rotation_matrix_SAMS);
		if (iResult == 0) // This happens if SIMA simulates more timesteps than SAMS
		{
			printf("Warning: SAMS connection is closed\n");
			(*ierr)++;
		}
		else if (iResult < 0)
		{
			printf("Receive failed from SAMS: %d\n", WSAGetLastError());
			*ierr = -29;
			return;
		}
		// On the first timestep, ice forces are absent.
		// The state-dependent forces Mx''+Kx are equal to the sea forces. F_coupled = F_sea
		if (send_to_SAMS(sams_tcp_socket, 0.0, F_coupled_local) == SOCKET_ERROR)
		{
			printf("Error sending TCP data to SAMS: %d\n", WSAGetLastError());
			closesocket(sams_tcp_socket);
			WSACleanup();
			*ierr = -28;
			return;
		}

		char itconfig_file_path[MCHEXT]; // Also used to get the SAMS simulation name

		// From the given folder, find the SAMS result file with the latest timestamp in the name
		{
			char SAMS_resultfolder_path[MCHEXT];
			char* string_end;
			char* SAMS_basename_start;
			char* SAMS_basename_end;
			int SAMS_basename_length;
			int SAMS_basename_offset;
			char SAMS_basename[MCHEXT];
			char SAMS_result_mask[MAX_PATH];
			WIN32_FIND_DATA found_data;
			HANDLE found_handle = NULL;
			char SAMS_basename_date_mask[MCHEXT];
			int latest_year = 0;
			int latest_month = 0;
			int latest_day = 0;
			int latest_time = 0;
			int SAMS_result_year;
			int SAMS_result_month;
			int SAMS_result_day;
			int SAMS_result_time;
			char SAMS_latest_result[MCHEXT];
			bool newer_result_found;
			// These string parameters are defined in the "External DLL force" dialog in SIMA
			strncpy(itconfig_file_path, chext[0] + 1, MCHEXT - 1);
			strncpy(SAMS_resultfolder_path, chext[2] + 1, MCHEXT - 1);
			// They are passed to this function with a leading space.
			// By the way, every parameter in sys.dat is recorded with a leading space
			// Also, the strings are given in contiguous fixed-width, space-padded blocks without zero-termination
			// Add zero-termination to the end of every string
			itconfig_file_path[MCHEXT - 1] = '\0';
			SAMS_resultfolder_path[MCHEXT - 1] = '\0';
			// Trim padding space
			string_end = SAMS_resultfolder_path + strlen(SAMS_resultfolder_path) - 1;
			while (string_end > SAMS_resultfolder_path && isspace((unsigned char)*string_end)) string_end--;
			string_end[1] = '\0'; // Write new null terminator character
			// The basename begins after the last backslash in the path
			SAMS_basename_start = strrchr(itconfig_file_path, '\\') + 1;
			// The basename ends before the file extension
			SAMS_basename_end = strstr(itconfig_file_path, ".itconfig");
			SAMS_basename_length = SAMS_basename_end - SAMS_basename_start;
			SAMS_basename_offset = SAMS_basename_start - itconfig_file_path;
			strncpy(SAMS_basename, itconfig_file_path + SAMS_basename_offset, SAMS_basename_length);
			SAMS_basename[SAMS_basename_length] = '\0';
			sprintf(SAMS_result_mask, "%s%s*.txt", SAMS_resultfolder_path, SAMS_basename);
			// Find all .txt files at the result folder matching the basename
			if ((found_handle = FindFirstFile(SAMS_result_mask, &found_data)) == INVALID_HANDLE_VALUE)
			{
				printf("Path not found when looking for SAMS result files: %s\n", SAMS_result_mask);
				*ierr = -1;
				return;
			}
			sprintf(SAMS_basename_date_mask, "%s_ %%2d%%2d%%4d_%%6d", SAMS_basename);
			do
			{
				newer_result_found = 0;
				iResult = sscanf
				(
					found_data.cFileName,
					SAMS_basename_date_mask,
					&SAMS_result_day,
					&SAMS_result_month,
					&SAMS_result_year,
					&SAMS_result_time
				);
				if (iResult != 4)
				{
					printf("Error parsing timestamps from SAMS result files\n");
					*ierr = 1;
					return;
				}
				if (SAMS_result_year > latest_year)
					newer_result_found = 1;
				else if (SAMS_result_year == latest_year)
				{
					if (SAMS_result_month > latest_month)
						newer_result_found = 1;
					else if (SAMS_result_month == latest_month)
					{
						if (SAMS_result_day > latest_day)
							newer_result_found = 1;
						else if (SAMS_result_day == latest_day)
							if (SAMS_result_time > latest_time)
								newer_result_found = 1;
					}

				}
				if (newer_result_found)
				{
					latest_year = SAMS_result_year;
					latest_month = SAMS_result_month;
					latest_day = SAMS_result_day;
					latest_time = SAMS_result_time;
					strcpy(SAMS_latest_result, found_data.cFileName);
				}
			} while (FindNextFile(found_handle, &found_data)); // Find the next file.
			FindClose(found_handle); // Always, Always, clean things up!
			sprintf(SAMS_resultfile_path, "%s%s", SAMS_resultfolder_path, SAMS_latest_result);
			printf("Reading SAMS result file: %s\n", SAMS_resultfile_path);
			SAMS_result_file = fopen(SAMS_resultfile_path, "r");
			if (!SAMS_result_file)
			{
				printf("Error opening SAMS result file\n");
				*ierr = -2;
				return;
			}
			if (SAMS_result_file == NULL)
			{
				printf("Timestep %d: Error: SAMS result file handle is NULL\n", step_nr);
				*ierr = -3;
				return;
			}
		};
		// Write the log headers
		{
			csv_writer = CsvWriter_new(gfexfo_result_file_name, ";", 0);
			for (int i = 0; i < nr_of_csv_ints; i++)
			{
				if (CsvWriter_writeField(csv_writer, csv_ints_header[i]))
				{
					printf("Error writing .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
					*ierr = -4;
					return;
				}
			}
			for (int i = 0; i < nr_of_csv_floats; i++)
			{
				if (CsvWriter_writeField(csv_writer, csv_floats_header[i]))
				{
					printf("Error writing to .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
					*ierr = -5;
					return;
				}
			}
			for (int i = 0; i < nr_of_csv_doubles; i++)
			{
				if (CsvWriter_writeField(csv_writer, csv_doubles_header[i]))
				{
					printf("Error writing .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
					*ierr = -6;
					return;
				}
			}
		}
		// Read the M and K matrices from sys.dat
		{
			FILE* sys_dat_file;
			char* sys_dat_filename = "sys-sima.dat";
			sys_dat_file = fopen(sys_dat_filename, "r");
			if (!sys_dat_file)
			{
				sys_dat_filename = "sys.dat";
				sys_dat_file = fopen(sys_dat_filename, "r");
				if (!sys_dat_file)
				{
					printf("Could not open either sys.dat or sys-sima.dat\n");
					*ierr = -7;
					return;
				}
			}
			while (getline(&line_buffer, &line_buffer_size, sys_dat_file) >= 0)
			{
				if (strcmp(line_buffer, " MASS COEFFICIENTS\n") == 0)
					break;
			}
			// Skip the separator line and headers
			for (int i = 0;i < 2;i++)
				getline(&line_buffer, &line_buffer_size, sys_dat_file);
			if (getline(&line_buffer, &line_buffer_size, sys_dat_file) != 114)
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
				*ierr = -8;
				return;
			}
			mass_matrix_SIMA[1][1] = mass_matrix_SIMA[0][0];
			mass_matrix_SIMA[2][2] = mass_matrix_SIMA[0][0];
			while (getline(&line_buffer, &line_buffer_size, sys_dat_file) >= 0)
			{
				if (strcmp(line_buffer, " STIFFNESS REFERENCE\n") == 0)
					break;
			}
			for (int i = 0;i < 2;i++)
				getline(&line_buffer, &line_buffer_size, sys_dat_file); // Skip the separator line and headers

			int line_length = getline(&line_buffer, &line_buffer_size, sys_dat_file);
			if (line_length < 98 || line_length > 104)
			printf("Warning reading stiffness reference from %s: unexpected line length\n", sys_dat_filename);
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
				*ierr = -9;
				return;
			}
			while (getline(&line_buffer, &line_buffer_size, sys_dat_file) >= 0)
			{
				if (strcmp(line_buffer, "'KMAT\n") == 0)
					break;
			}
			for (int i = 0;i < 6;i++)
			{
				if (getline(&line_buffer, &line_buffer_size, sys_dat_file) != 86)
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
					*ierr = -10;
					return;
				}
			}
			fclose(sys_dat_file);
		};
		// Check that the simulation is set up similarly in SAMS
		{
			// TODO this search routine is done more than once, rewrite as function
			FILE* itconfig_file = fopen(itconfig_file_path, "r");
			if (!itconfig_file)
			{
				printf("Error opening:%s\n", itconfig_file_path);
				*ierr = -11;
				return;
			}
			while (getline(&line_buffer, &line_buffer_size, itconfig_file) >= 0)
			{
				if (strstr(line_buffer, "\"maxSimulationTime\""))
					break;
			}
			double SAMS_maxSimulationTime;
			if (sscanf(line_buffer, " \" maxSimulationTime \" : %lf ", &SAMS_maxSimulationTime) != 1)
			{
				printf("Error parsing maxSimulationTime from:%s\n", itconfig_file_path);
				*ierr = -12;
				return;
			}
			if (SAMS_maxSimulationTime != nr_of_steps  * dt)
			{
				printf("Error: Different simulation times declared in SAMS and SIMA.\n");
				printf("SAMS maxSimulationTime : %.1f SIMA : %.1f dt: %.1f\n", SAMS_maxSimulationTime, nr_of_steps * dt, dt);
				*ierr = -30;
				return;
			}
			// TODO also check that the frequency is the same, in .itconfig it is given as f=1/dt
			fclose(itconfig_file);
		};
		double radius_of_gyration_SAMS[3];
		// Read the structure mass coefficients from SAMS
		{
			char structure_file_path[MCHEXT];
			strncpy(structure_file_path, chext[1] + 1, MCHEXT - 1);
			structure_file_path[MCHEXT - 1] = '\0';
			FILE* structure_file = fopen(structure_file_path, "r");
			if (!structure_file)
			{
				printf("Error opening:%s\n", structure_file_path);
				*ierr = -13;
				return;
			}
			while (getline(&line_buffer, &line_buffer_size, structure_file) >= 0)
			{
				if (strstr(line_buffer, "\"structureMass\""))
					break;
			}
			if (sscanf(line_buffer, " \" structureMass \" : %lf ", &mass_matrix_SAMS[0][0]) != 1)
			{
				printf("Error parsing structure mass from:%s\n", structure_file_path);
				*ierr = -14;
				return;
			}
			mass_matrix_SAMS[1][1] = mass_matrix_SAMS[0][0];
			mass_matrix_SAMS[2][2] = mass_matrix_SAMS[0][0];
			// TODO Also check that the frequency is same as SIMA. In .structure.txt, it is given as f=1/dt
			// Do not know whether the mass or the radii come first, so search from the top
			fclose(structure_file);
			structure_file = fopen(structure_file_path, "r");
			while (getline(&line_buffer, &line_buffer_size, structure_file) >= 0)
			{
				if (strstr(line_buffer, "\"structureRadiusOfGyration\""))
					break;
			}
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
				*ierr = -15;
				return;
			}
			fclose(structure_file);
		}
		// Calculate the mass moments of inertia from the radii of gyration
		for (int i = 0;i < 3;i++)
			mass_matrix_SAMS[3 + i][3 + i] = pow(radius_of_gyration_SAMS[i], 2) * mass_matrix_SAMS[i][i];
		// Check the mass matrix for disrepancies between SAMS and SIMA
		double mass_matrix_relative_difference[6][6];
		for (int i = 0;i < 6;i++)
		{
			for (int j = 0;j < 6;j++)
			{
				if (mass_matrix_SAMS[i][j] && mass_matrix_SIMA[i][j])
				{
					mass_matrix_relative_difference[i][j] = (mass_matrix_SIMA[i][j] * 1000 - mass_matrix_SAMS[i][j]) / mass_matrix_SIMA[i][j];
					if (mass_matrix_relative_difference[i][j] > 0.01)
					{
						printf
						(
							"Warning: At position [%d][%d], the mass matrix has %lf in SIMA and %lf in SAMS, relative difference %le\n",
							i,
							j,
							mass_matrix_SIMA[i][j] * 1000,
							mass_matrix_SAMS[i][j],
							mass_matrix_relative_difference[i][j]
						);
						(*ierr)++;
					}
				}
			}
		}
	};

	/*fclose(SAMS_result_file);
	SAMS_result_file = fopen(SAMS_resultfile_path, "r");
	int line_counter = 0;
	while (true)
	{
		int line_length = getline(&line_buffer, &line_buffer_size, SAMS_result_file);
		if (line_length < 0)
			break;
		line_counter++;
	}
	printf("Step %d: %d lines read\n", step_nr, line_counter);
	printf("Last line: %s\n",line_buffer);*/

	//int lines_to_read;
	//if (step_nr < 3) // Headers not written yet, avoid EOF error
	//	lines_to_read = 0;
	//else if (step_nr == 3) // Headers written , catch up
	//	lines_to_read = 7;
	//else // Keep catching up
	//	lines_to_read = 1;
	//for (int i = 0;i < lines_to_read;i++)
	//{
	//	if (SAMS_result_file == NULL)
	//	{
	//		printf("Timestep %d: Error: SAMS result file handle is NULL\n",step_nr);
	//		*ierr = -17;
	//		return;
	//	}
	//	int line_length = getline(&line_buffer, &line_buffer_size, SAMS_result_file);
	//	
	//	if (line_length < 0)
	//	{
	//		printf("Timestep %d: Error getting line %d from SAMS result file\n",step_nr,i);
	//		*ierr = -18;
	//		return;
	//	}
	//	printf("Step %d, line %d: %s\n", step_nr, i, line_buffer);
	//}
	//if (lines_to_read)
	//{
	//	int current_index = 0;
	//	int offset = 0;

	//	for (int i = 0;i < 47;i++) // 47 columns in the SAMS text output file
	//	{
	//		if (sscanf(line_buffer + current_index, " %lf %n", &SAMS_txt_output[i], &offset) != 1)
	//		{
	//			// TODO show the actual file paths in these messages
	//			printf("Error parsing column %d from SAMS result file\n",i);
	//			//*ierr = -19;
	//			//return;
	//		}
	//		current_index += offset;
	//	}
	//	// Rename some SAMS text output values for clarity
	//	double SAMS_txt_time = SAMS_txt_output[0];
	//	if (SAMS_txt_time - time >= dt)
	//	{
	//		printf("Time mismatch: SIMA %f, SAMS log %f\n", time, SAMS_txt_time);
	//		*ierr = -20;
	//		return;
	//	}
	//	for (int i = 0;i < 3;i++)
	//	{
	//		// Linear forces are composed of breaking and rubble forces
	//		SAMS_ice_force[i] = SAMS_txt_output[20 + i] + SAMS_txt_output[23 + i];
	//		// Rotational forces are given as total
	//		SAMS_ice_force[3 + i] = SAMS_txt_output[26 + i];
	//	}
	//}

	// Last timestep, SAMS won't send anything. Close the connection and logs
	if (step_nr == nr_of_steps && substep_nr == nr_of_substeps)
	{
		iResult = WSACleanup();
		if (iResult == SOCKET_ERROR)
		{
			printf("Error closing connection: %d\n", WSAGetLastError());
			*ierr = -26;
			return;
		}
		CsvWriter_destroy(csv_writer);
		printf("Disconnected from SAMS, results saved in %s\n", gfexfo_result_file_name);
		fclose(SAMS_result_file); // TODO clean up this file also in case of ungraceful exit
		return;
	}

	// Receive and process data from SAMS
	{
		if (sams_tcp_socket == INVALID_SOCKET)
		{
			printf("Invalid SAMS TCP socket\n");
			*ierr = -16;
			return;
		}
		
		iResult = receive_from_SAMS(sams_tcp_socket, &SAMS_time, displacement_SAMS, velocity_SAMS_t0, rotation_matrix_SAMS);
		if (iResult == 0) // This happens if SIMA simulates more timesteps than SAMS
		{
			printf("Warning: SAMS connection is closed\n");
			(*ierr)++;
		}
		else if (iResult < 0)
		{
			printf("Receive failed from SAMS: %d\n", WSAGetLastError());
			*ierr = -21;
			return;
		}
		// Additional analysis on data received from SAMS, a bit redundant
		for (int i = 0;i < 6;i++)
		{
			// Calculate accelerations as backward finite differences of velocity.
			acceleration_SAMS[i] =
				(
					-velocity_SAMS_tm1[i]
					+ velocity_SAMS_t0[i]
					) / dt;
			// Calculate inertial forces
			inertial_force_SAMS[i] = acceleration_SAMS[i] * mass_matrix_SAMS[i][i];
			// Calculate hydrostatic forces, see supervision presentation of 27.04.2022
			hydrostatic_force_SAMS[i] = (stiffness_reference_SIMA[i] - displacement_SAMS[i]) * stiffness_matrix_SIMA[i][i];
		}
	}
	// Send data to SAMS
	if (send_to_SAMS(sams_tcp_socket, time, F_coupled_local) == SOCKET_ERROR)
	{
		printf("Error sending TCP data to SAMS: %d\n", WSAGetLastError());
		closesocket(sams_tcp_socket);
		WSACleanup();
		*ierr = -22;
		return;
	}

	// Log the variables for this timestep
	{
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
			SAMS_time, // SAMS_Time
			displacement_SAMS[0],
			displacement_SAMS[1],
			displacement_SAMS[2],
			displacement_SAMS[3],
			displacement_SAMS[4],
			displacement_SAMS[5],
			acceleration_SAMS[0],
			acceleration_SAMS[1],
			acceleration_SAMS[2],
			acceleration_SAMS[3],
			acceleration_SAMS[4],
			acceleration_SAMS[5],
			inertial_force_SAMS[0],
			inertial_force_SAMS[1],
			inertial_force_SAMS[2],
			inertial_force_SAMS[3],
			inertial_force_SAMS[4],
			inertial_force_SAMS[5],
			hydrostatic_force_SAMS[0],
			hydrostatic_force_SAMS[1],
			hydrostatic_force_SAMS[2],
			hydrostatic_force_SAMS[3],
			hydrostatic_force_SAMS[4],
			hydrostatic_force_SAMS[5],
			rotation_matrix_SIMA[0][0],
			rotation_matrix_SIMA[0][1],
			rotation_matrix_SIMA[0][2],
			rotation_matrix_SIMA[1][0],
			rotation_matrix_SIMA[1][1],
			rotation_matrix_SIMA[1][2],
			rotation_matrix_SIMA[2][0],
			rotation_matrix_SIMA[2][1],
			rotation_matrix_SIMA[2][2],
			gamma_float.double_precision[0],
			gamma_float.double_precision[1],
			gamma_float.double_precision[2],
			gamma_float.double_precision[3],
			rotation_matrix_SAMS[0][0],
			rotation_matrix_SAMS[0][1],
			rotation_matrix_SAMS[0][2],
			rotation_matrix_SAMS[1][0],
			rotation_matrix_SAMS[1][1],
			rotation_matrix_SAMS[1][2],
			rotation_matrix_SAMS[2][0],
			rotation_matrix_SAMS[2][1],
			rotation_matrix_SAMS[2][2],
			SAMS_ice_force[0],
			SAMS_ice_force[1],
			SAMS_ice_force[2],
			SAMS_ice_force[3],
			SAMS_ice_force[4],
			SAMS_ice_force[5]
		};
		CsvWriter_nextRow(csv_writer);
		char log_buffer[50];
		for (int i = 0; i < nr_of_csv_ints; i++)
		{
			sprintf(log_buffer, "%d", csv_ints[i]);
			if (CsvWriter_writeField(csv_writer, log_buffer))
			{
				printf("Error writing .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
				*ierr = -23;
				return;
			}
		}
		for (int i = 0; i < nr_of_csv_floats; i++)
		{
			sprintf(log_buffer, "%f", csv_floats[i]);
			if (CsvWriter_writeField(csv_writer, log_buffer))
			{
				printf("Error writing .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
				*ierr = -24;
				return;
			}
		}
		for (int i = 0; i < nr_of_csv_doubles; i++)
		{
			sprintf(log_buffer, "%f", csv_doubles[i]);
			if (CsvWriter_writeField(csv_writer, log_buffer))
			{
				printf("Error writing .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
				*ierr = -25;
				return;
			}
		}
	}

	// Calculations done, prepare for the next timestep
	for (int i = 0; i < 6; i++)
	{
		displacement_SIMA_tm2[i] = displacement_SIMA_tm1[i];
		displacement_SIMA_tm1[i] = displacement_SIMA_t0.double_precision[i];
		velocity_SAMS_tm1[i] = velocity_SAMS_t0[i];
	}

	free(line_buffer); // Not clear why I have to free this and not everything else
	//fflush(NULL); // supposed to write all buffers to file, but actually doesn't
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
		return INVALID_SOCKET;
}

int send_to_SAMS(SOCKET sams_tcp_socket, double time, double central_forces[6])
{
	// Only central forces implemented here, any point forces acting elsewhere on the body will not be sent
	// Encode doubles into a bytestring for SAMS
	union
	{
		double double_precision[1 + 6]; // time and central forces in 6DOF
		char character_array[(1 + 6) * sizeof(double)];
	} GFEXFO_to_SAMS;
	// SIMO receives forces through stor[] in kN and MN*m
	// SAMS receives forces through TCP in N and N*m
	// TODO leave unit consistency to the SIMA user, there are SIMA options, check them if possible
	GFEXFO_to_SAMS.double_precision[0] = time; // [s]
	for (int i = 0;i < 6;i++)
		GFEXFO_to_SAMS.double_precision[1 + i] = central_forces[i];
	return send(sams_tcp_socket, GFEXFO_to_SAMS.character_array, sizeof(GFEXFO_to_SAMS.double_precision), 0);
}

int receive_from_SAMS(SOCKET sams_tcp_socket, double* SAMS_time, double displacement_SAMS[6], double velocity_SAMS_t0[6], double rotation_matrix_SAMS[3][3])
{
	// Decode a bytestring from SAMS into doubles
	union
	{
		double double_precision[32];
		char character_array[32 * sizeof(double)];
	} SAMS_to_GFEXFO;
	// SAMS writes output before it reads input (text based used control features of SAMS) 
	int iResult = recv(sams_tcp_socket, SAMS_to_GFEXFO.character_array, sizeof(SAMS_to_GFEXFO.double_precision), 0);
	if (iResult <= 0)
		return iResult;
	*SAMS_time = SAMS_to_GFEXFO.double_precision[0];
	// Rename the axis components and the angle for clarity
	double e1 = SAMS_to_GFEXFO.double_precision[4];
	double e2 = SAMS_to_GFEXFO.double_precision[5];
	double e3 = SAMS_to_GFEXFO.double_precision[6];
	double theta_aa = SAMS_to_GFEXFO.double_precision[7]; // Theta in axis-angle representation
	// Calculations based on data received from SAMS
	rotation_matrix_SAMS[0][0] = (1 - cos(theta_aa)) * e1 * e1 + cos(theta_aa);
	rotation_matrix_SAMS[0][1] = (1 - cos(theta_aa)) * e1 * e2 - e3 * sin(theta_aa);
	rotation_matrix_SAMS[0][2] = (1 - cos(theta_aa)) * e1 * e3 + e2 * sin(theta_aa);
	rotation_matrix_SAMS[1][0] = (1 - cos(theta_aa)) * e2 * e1 + e3 * sin(theta_aa);
	rotation_matrix_SAMS[1][1] = (1 - cos(theta_aa)) * e2 * e2 + cos(theta_aa);
	rotation_matrix_SAMS[1][2] = (1 - cos(theta_aa)) * e2 * e3 - e1 * sin(theta_aa);
	rotation_matrix_SAMS[2][0] = (1 - cos(theta_aa)) * e3 * e1 - e2 * sin(theta_aa);
	rotation_matrix_SAMS[2][1] = (1 - cos(theta_aa)) * e3 * e2 + e1 * sin(theta_aa);
	rotation_matrix_SAMS[2][2] = (1 - cos(theta_aa)) * e3 * e3 + cos(theta_aa);
	// TODO express the difference between SAMS and SIMA rotation matrices, define a maximum value (1%?)
	// Convert the Euler axis-angle representation to Tait-Bryan chained Z-Y-X rotations
	// The calculations were analytically derived from the rotation matrix definition in SIMO manual, app.C
	double theta_zyx = asin(-rotation_matrix_SAMS[2][0]); // Theta in Z-Y-X representation
	double phi = asin(rotation_matrix_SAMS[2][1] / cos(theta_zyx));
	double psi = asin(rotation_matrix_SAMS[1][0] / cos(theta_zyx));
	displacement_SAMS[3] = phi;
	displacement_SAMS[4] = theta_zyx;
	displacement_SAMS[5] = psi;
	for (int i = 0;i < 3;i++)
		displacement_SAMS[i] = SAMS_to_GFEXFO.double_precision[1 + i];
	for (int i = 0;i < 6;i++)
		velocity_SAMS_t0[i] = SAMS_to_GFEXFO.double_precision[8 + i];
	return iResult;
}