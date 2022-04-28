#define _USE_MATH_DEFINES
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
	static double displacement_SIMA_tm1[6]; // Timestep t = -1.
	static double displacement_SIMA_tm2[6]; // Timestep t = -2
	static double velocity_SAMS_tm1[6]; // Timestep t = -1
	static double mass_matrix_SIMA[6][6]; // [t=kg*10^3]
	static double mass_matrix_SAMS[6][6]; // [kg]
	static double stiffness_matrix_SIMA[6][6]; // Hydrostatic restoring force is modeled as spring stiffness
	static double stiffness_reference_SIMA[6]; // Hydrostatic equilibrium position
	static char SAMS_resultfile_path[MCHEXT];

	// Assign constant forces
	if (time > 5.)
		stor[0] = 10.; // [kN] surge
	else
		stor[0] = 0.;
	stor[1] = 0.; // [kN] sway
	stor[2] = 0.; // [kN] heave
	stor[3] = 0.; // [MN*m] roll
	stor[4] = 0.; // [MN*m] pitch
	if (time >= 0.5 && time < 4.5)
		stor[5] = 0.1; // [MN*m] yaw
	else
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
		"displacement_SAMS_0",
		"displacement_SAMS_1",
		"displacement_SAMS_2",
		"displacement_SAMS_3",
		"displacement_SAMS_4",
		"displacement_SAMS_5",
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
	};
	int nr_of_csv_ints = (int)sizeof(csv_ints_header) / sizeof(csv_ints_header[0]);
	int nr_of_csv_floats = (int)sizeof(csv_floats_header) / sizeof(csv_floats_header[0]);
	int nr_of_csv_doubles = (int)sizeof(csv_doubles_header) / sizeof(csv_doubles_header[0]);
	char* result_file_name = "gfexfo_sams.csv";
	
	// For error codes
	int iResult;


	
	run_counter++;

	double velocity_SIMA_t0[6];
	double velocity_SAMS_t0[6];
	double acceleration_SIMA[6];
	double displacement_SAMS[6];
	double acceleration_SAMS[6];
	double inertial_force_SIMA[6]; // [kN, MN*m]
	double inertial_force_SAMS[6]; // [N]
	double hydrostatic_force_SIMA[6]; // [kN, MN*m]
	double hydrostatic_force_SAMS[6]; // [N]
	double radius_of_gyration_SAMS[3];
	double mass_matrix_difference[6][6];
	double rotation_matrix_SIMA[3][3];
	double rotation_matrix_SAMS[3][3];
	double rotation_axis_SAMS[3];
	double rotation_angle_SAMS;
	double rotation_velocity_negative_jump;
	double rotation_velocity_positive_jump;
	double rotation_acceleration_negative_jump;
	double rotation_acceleration_positive_jump;

	// Back-initialize the SIMA displacements to avoid initial jerk
	if (run_counter == 1)
	{
		for (int i = 0;i < 6;i++)
		{
			displacement_SIMA_tm1[i] = displacement_SIMA_t0.double_precision[i];
			displacement_SIMA_tm2[i] = displacement_SIMA_t0.double_precision[i];
		}
	}

	// First run, open the SAMS connection and write the log headers
	if (step_nr == 1 && substep_nr == 1 && iteration_nr == 0)
	{

		char itconfig_file_path[MCHEXT];
		// These string parameters are defined in the "External DLL force" dialog in SIMA
		// They are passed to this function with a leading space and without a zero termination
		strncpy(itconfig_file_path, chext[0] + 1, MCHEXT - 1);
		itconfig_file_path[MCHEXT - 1] = '\0';
		char structure_file_path[MCHEXT];
		strncpy(structure_file_path, chext[1] + 1, MCHEXT - 1);
		structure_file_path[MCHEXT - 1] = '\0';
		char SAMS_resultfolder_path[MCHEXT];
		strncpy(SAMS_resultfolder_path, chext[2] + 1, MCHEXT - 1);
		// Trim trailing space
		char* string_end = SAMS_resultfolder_path + strlen(SAMS_resultfolder_path) - 1;
		while (string_end > SAMS_resultfolder_path && isspace((unsigned char)*string_end)) string_end--;
		string_end[1] = '\0'; // Write new null terminator character
		// The basename begins after the last backslash in the path
		char* SAMS_basename_start = strrchr(itconfig_file_path, '\\') + 1;
		// The basename ends before the file extension
		char* SAMS_basename_end = strstr(itconfig_file_path, ".itconfig");
		int SAMS_basename_length = SAMS_basename_end - SAMS_basename_start;
		int SAMS_basename_offset = SAMS_basename_start - itconfig_file_path;
		char SAMS_basename[MCHEXT];
		strncpy(SAMS_basename, itconfig_file_path + SAMS_basename_offset, SAMS_basename_length);
		SAMS_basename[SAMS_basename_length] = '\0';
		char SAMS_result_mask[MAX_PATH];
		sprintf(SAMS_result_mask, "%s%s*.txt", SAMS_resultfolder_path, SAMS_basename);
		// Find all .txt files at the result folder matching the basename
		WIN32_FIND_DATA found_data;
		HANDLE found_handle = NULL;
		if ((found_handle = FindFirstFile(SAMS_result_mask, &found_data)) == INVALID_HANDLE_VALUE)
		{
			printf("Path not found: %s\n", SAMS_result_mask);
			*ierr=1;
			return;
		}
		char SAMS_basename_date_mask[MCHEXT];
		sprintf(SAMS_basename_date_mask, "%s_ %%2d%%2d%%4d_%%6d", SAMS_basename);
		int latest_year=0;
		int latest_month=0;
		int latest_day=0;
		int latest_time=0;
		int SAMS_result_year;
		int SAMS_result_month;
		int SAMS_result_day;
		int SAMS_result_time;
		char SAMS_latest_result[MCHEXT];
		bool newer_result_found;
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
		} while (FindNextFile(found_handle, &found_data)); //Find the next file.
		FindClose(found_handle); //Always, Always, clean things up!
		
		sprintf(SAMS_resultfile_path, "%s%s", SAMS_resultfolder_path, SAMS_latest_result);
		printf("Latest: %s\n", SAMS_resultfile_path);
		FILE* SAMS_result_file = fopen(SAMS_resultfile_path,"r");
		if (!SAMS_result_file)
		{
			printf("Error opening SAMS result file\n");
			*ierr = 1;
			return;
		}
		fclose(SAMS_result_file);
		printf("SAMS result file succesfully opened and closed\n");

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
		// TODO check that SAMS simulates one timestep less than SIMA
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
		// Check the mass matrix for disrepancies between SAMS and SIMA
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
	for (int i = 0;i < 3;i++)
	{
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
	for (int i = 0;i < 6;i++)
	{
		// Calculate linear forces
		inertial_force_SIMA[i] = acceleration_SIMA[i] * mass_matrix_SIMA[i][i];
		hydrostatic_force_SIMA[i] =
			(stiffness_reference_SIMA[i] - displacement_SIMA_t0.double_precision[i])
			* stiffness_matrix_SIMA[i][i];
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

		// SAMS writes output before it reads input (text based used control features of SAMS) 
		// Receive kinematic and dynamic information from SAMS
		iResult = recv(sams_tcp_socket, SAMS_to_GFEXFO.character_array, sizeof(SAMS_to_GFEXFO.double_precision), 0);
		// TODO get additional data from SAMS log
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
		{
			// Rename some received variables for clarity
			displacement_SAMS[i] = SAMS_to_GFEXFO.double_precision[1 + i];
		}
		for (int i = 0;i < 6;i++)
		{
			velocity_SAMS_t0[i] = SAMS_to_GFEXFO.double_precision[8 + i];
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

		// SIMO receives forces through stor[] in kN and kN*m
		// SAMS receives forces through TCP in N and N*m
		GFEXFO_to_SAMS.double_precision[0] = (double)time; // [s]
		GFEXFO_to_SAMS.double_precision[1] = stor[0] * pow(10, 3); // 3 to get the same resultant acceleration
		GFEXFO_to_SAMS.double_precision[2] = stor[1] * pow(10, 3);
		GFEXFO_to_SAMS.double_precision[3] = stor[2] * pow(10, 3);
		GFEXFO_to_SAMS.double_precision[4] = stor[3] * pow(10, 6); // 6 to get the same resultant acceleration. Why though?
		GFEXFO_to_SAMS.double_precision[5] = stor[4] * pow(10, 6);
		GFEXFO_to_SAMS.double_precision[6] = stor[5] * pow(10, 6);
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
			displacement_SAMS[0],
			displacement_SAMS[1],
			displacement_SAMS[2],
			displacement_SAMS[3],
			displacement_SAMS[4],
			displacement_SAMS[5],
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
			rotation_matrix_SAMS[2][2]
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
