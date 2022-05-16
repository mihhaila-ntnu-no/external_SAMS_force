#define _USE_MATH_DEFINES
#include <math.h>
#include <stdbool.h>
#include <direct.h>
#include "simo_extfunc.h"
#include "csvwriter.h"
#include "getline.h"
#include "errno.h"

extern int errno;

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
	// Lowest used value: -34
	*ierr = 0;

	// TODO try simulating with more than 1 substep, is dt referring to the length of one step or substep?
	// Function argument aliases
	int body_nr = iinfo[2]-1; // bodies are indexed starting with 1
	int step_nr = iinfo[5];
	int nr_of_steps = iinfo[6];
	int substep_nr = iinfo[7];
	int nr_of_substeps = iinfo[8];
	int substep_nr_overall = substep_nr + (step_nr - 1) * nr_of_substeps;
	int nr_of_substeps_overall = nr_of_steps * nr_of_substeps;
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
	static double mass_matrix_SIMA[6][6]; // [kg*10^3, kg*m^2*10^3]
	static double mass_matrix_SAMS[6][6]; // [kg]
	static char SAMS_resultfile_path[MCHEXT];
	static int nr_of_csv_ints; // For logging
	static int nr_of_csv_floats;
	static int nr_of_csv_doubles;
	static bool called_by_RIFLEX = true; // First assumption
	static bool called_by_SIMO = false;
	static double F_ice_SIMA_local[6];
	// Hydrostatic restoring force is modeled as spring stiffness
	static double stiffness_matrix_SIMA[6][6]; // [kN/m, kN*m]
	static double stiffness_reference_SIMA[6]; // [m, rad] Hydrostatic equilibrium position

	// If this is not the first iteration of the timestep, return the forces from the first one
	if (iinfo[11] > 0) // iinfo[11] = iteration nr, starting from 0
	{
		for (int i = 0;i < 3;i++)
		{
			stor[i] = F_ice_SIMA_local[i];
			stor[3 + i] = F_ice_SIMA_local[3 + i];
			stor[6 + i] = 0.; // "internal parameter", don't know what for
		}
		*ierr = 0;
		return;
	}

	// Working variables
	double velocity_SIMA_t0[6];
	double acceleration_SIMA[6];
	double displacement_SAMS[6];
	double acceleration_SAMS[6];
	double inertial_force_SIMA[6]; // [kN, MN*m]
	double inertial_force_SAMS[6]; // [N]
	double hydrostatic_force_SIMA[6]; // [kN, MN*m]
	double hydrostatic_force_SAMS[6]; // [N]
	double rotation_matrix_SIMA[3][3]; // From global to local
	double rotation_matrix_SAMS[3][3]; // So X_global = matrix x X_local
	double F_ice_global[6];
	double F_sea_global[6]; // Sea forces on structure, as opposed to ice forces
	double F_sea_SAMS_local[6];
	double F_coupled_global[6]; // [kN, MN*m] Forces from both sea and ice, calculated by SIMA
	double timestep_start_time;
	double timestep_end_time;
	double SAMS_TCP_time;
	double SAMS_txt_time;
	
	int iResult; // For error codes
	char* line_buffer = NULL; // for getline()
	size_t line_buffer_size = 0; // for getline()

	char* gfexfo_result_file_name = "gfexfo_sams.csv";
	char itconfig_file_path[MCHEXT]; // Also used to get the SAMS simulation name
	
	// Simulation setup. Read and assume information.
	if (run_counter == 1)
	{
		// Back-initialize the SIMA displacements to avoid calculating an unphysical initial jerk
		for (int i = 0;i < 6;i++)
		{
			displacement_SIMA_tm1[i] = displacement_SIMA_t0.double_precision[i];
			displacement_SIMA_tm2[i] = displacement_SIMA_t0.double_precision[i];
		}
		// Read the M and K matrices from sys.dat, find out which SIMA module is calling this function
		{
			FILE* sys_dat_file;
			char* system_description_filename = "sys-sima.dat";
			sys_dat_file = fopen(system_description_filename, "r");
			if (!sys_dat_file)
			{
				// Second assumption: this function was called by SIMO
				called_by_RIFLEX = false;
				called_by_SIMO = true;
				system_description_filename = "sys.dat";
				sys_dat_file = fopen(system_description_filename, "r");
			}
			if (!sys_dat_file) // This function was called by neither RIFLEX nor SIMO
			{
				printf("Could not find system description file.\n");
				printf("A RIFLEX task should name it sys-sima.dat; a SIMO task should name it sys.dat\n");
				*ierr = -7;
				return;
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
				printf("Warning reading stiffness reference from %s: unexpected line length\n", system_description_filename);
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
		strncpy(itconfig_file_path, chext[0] + 1, MCHEXT - 1);
		itconfig_file_path[MCHEXT - 1] = '\0';
		FILE* itconfig_file = fopen(itconfig_file_path, "r");
		if (!itconfig_file)
		{
			printf("Error opening .itconfig file:%s\n", itconfig_file_path);
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
		// TODO check that the simulation time >= 19 s. Shorter simulations seem to cause SAMS to expect an extra TCP message at the end
		if (SAMS_maxSimulationTime != nr_of_steps * dt)
		{
			printf("Error: Different simulation times declared in SAMS and SIMA.\n");
			printf("SAMS maxSimulationTime : %.1f SIMA : %.1f dt: %.1f\n", SAMS_maxSimulationTime, nr_of_steps * dt, dt);
			*ierr = -30;
			return;
		}
		// TODO also check that the frequency is the same, in .itconfig it is given as f=1/dt
		// Here, dt should be the length on one substep, not step.
		// Otherwise, SAMS would be called only on the first substep of every step
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
				printf("Error opening .structure.txt file:%s\n", structure_file_path);
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
		// TODO add units to the warning text
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

		// TODO check that the hydrostatic properties are {} empty in the .structure.txt file

		// Connect to SAMS
		sams_tcp_socket = connect_to_SAMS();
		if (sams_tcp_socket == INVALID_SOCKET)
		{
			printf("Unable to connect to SAMS\n");
			WSACleanup();
			*ierr = -27;
			return;
		}

		// Find the SAMS result .txt file with the latest timestamp in the name
		{
			// This file will exist as soon as the SAMS simulation is launched, even if no information has been exchanged yet
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
			char SAMS_resultfile_name_latest[MCHEXT];
			bool newer_result_found;
			// Three such string parameters are defined in the "External DLL force" dialog in SIMA
			strncpy(SAMS_resultfolder_path, chext[2] + 1, MCHEXT - 1);
			// They are passed to this function with a leading space.
			// By the way, every parameter in sys.dat is recorded with a leading space
			// Also, the strings are given in contiguous fixed-width, space-padded blocks without zero-termination
			// Add zero-termination to the end of every string
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
					strcpy(SAMS_resultfile_name_latest, found_data.cFileName);
				}
			} while (FindNextFile(found_handle, &found_data)); // Find the next file.
			FindClose(found_handle); // Always, Always, clean things up!
			sprintf(SAMS_resultfile_path, "%s%s", SAMS_resultfolder_path, SAMS_resultfile_name_latest);
		}

		// Assume zero ice forces for the first timestep
		for (int i = 0;i < 6;i++)
			F_ice_global[i] = 0;

	};

	// Calculate kinematic state, state-dependent forces, global->body transformation matrix
	{
		if (called_by_RIFLEX)
		{
			timestep_start_time = rinfo[0] - dt;
			timestep_end_time = rinfo[0];
		}
		else if (called_by_SIMO)
		{
			timestep_start_time = rinfo[0];
			timestep_end_time = rinfo[0] + dt;
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
			// TODO these are rotational joint velocities, each of them has its own coordinate frame
			// This method of finding the accelerations neglects the fact that the coordinate frames change over dt
			// Rather, the joint velocities should be converted to a angular velocity vector in the global frame
			// Then, will its derivative be equal to an acceleration vector? Maybe, check literature
			// TODO a Google search even suggests that the acceleration vector could be found by differentiating the rotation matrix
			// With small angles, the rotational joint velocities are very close to the global-frame angular velocities
			// All this shown on paper notes, page 9
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
		// TODO units?
		for (int i = 0;i < 6;i++)
		{
			inertial_force_SIMA[i] = acceleration_SIMA[i] * mass_matrix_SIMA[i][i];
			hydrostatic_force_SIMA[i] =
				(stiffness_reference_SIMA[i] - displacement_SIMA_t0.double_precision[i])
				* stiffness_matrix_SIMA[i][i];
			F_coupled_global[i] = inertial_force_SIMA[i] + hydrostatic_force_SIMA[i];
		}
	};

	// Calculate sea forces
	for (int i = 0;i < 3;i++)
	{
		F_sea_global[i] = F_coupled_global[i] - F_ice_global[i] * pow(10, -3); // N -> kN
		F_sea_global[3+i] = F_coupled_global[3+i] - F_ice_global[3+i] * pow(10, -6); // N*m -> MN*m
	}
	
	// Receive TCP data from SAMS
	{
		if (sams_tcp_socket == INVALID_SOCKET)
		{
			printf("Invalid SAMS TCP socket\n");
			*ierr = -16;
			return;
		}
		iResult = receive_from_SAMS(sams_tcp_socket, &SAMS_TCP_time, displacement_SAMS, rotation_matrix_SAMS);
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
		if ((SAMS_TCP_time - timestep_start_time) >= dt / 10)
		{
			printf("Timestep start time mismatch: SIMA %f, SAMS TCP %f\n", timestep_start_time, SAMS_TCP_time);
			*ierr = -34;
			return;
		}
		// TODO express the difference between SAMS and SIMA rotation matrices, define a maximum value (1%?)
	}

	// Transform the sea forces from global coordinate system to SAMS local
	for (int i = 0;i < 3;i++)
	{
		F_sea_SAMS_local[i] = 0;
		F_sea_SAMS_local[3 + i] = 0;
		// Send zero forces until rhythm and coordinate transform are sorted
		//for (int j = 0;j < 3;j++)
		//{
		//	F_sea_SAMS_local[i] += rotation_matrix_SAMS[j][i] * F_sea_global[j] * pow(10, 3); // kN -> N
		//	F_sea_SAMS_local[3 + i] += rotation_matrix_SAMS[j][i] * F_sea_global[3 + j] * pow(10, 6); // MN*m -> N*m
		//}
	}

	// Send sea forces to SAMS
	if (send_to_SAMS(sams_tcp_socket, timestep_start_time, F_sea_SAMS_local) == SOCKET_ERROR)
	{
		printf("Error sending TCP data to SAMS: %d\n", WSAGetLastError());
		closesocket(sams_tcp_socket);
		WSACleanup();
		*ierr = -22;
		return;
	}

	// Read ice forces from the SAMS result .txt file
	{
		// For SIMO simulations, the i=0, t=0.0 line will not be written. Skip reading it
		iResult = read_from_SAMS(SAMS_resultfile_path, substep_nr_overall, nr_of_substeps_overall, &SAMS_txt_time, F_ice_global);
		if (iResult < 0)
		{
			*ierr = iResult;
			return;
		}
		if ((SAMS_txt_time - timestep_end_time) >= dt / 10)
		{
			printf("Timestep end time mismatch: SIMA %f, SAMS log %f\n", timestep_end_time, SAMS_txt_time);
			*ierr = -20;
			return;
		}
	}

	// Save the ice forces to SIMA's memory for the next timestep
	for (int i = 0;i < 3;i++)
	{
		// F_ice_SIMA_local holds the values from the previous timestep
		F_ice_SIMA_local[i] = 0; // [kN] surge, sway, heave
		F_ice_SIMA_local[3 + i] = 0; // [MN*m] roll, pitch, yaw
		// Transform the ice forces from global coordinate system to SIMA local
		for (int j = 0;j < 3;j++)
		{
			F_ice_SIMA_local[i] += rotation_matrix_SIMA[j][i] * F_ice_global[j] * pow(10, -3); // N -> kN
			F_ice_SIMA_local[3 + i] += rotation_matrix_SIMA[j][i] * F_ice_global[3 + j] * pow(10, -6); // N*m -> MN*m
		}
		// Save ice forces to SIMA memory
		stor[i] = F_ice_SIMA_local[i];
		stor[3 + i] = F_ice_SIMA_local[3 + i];
		stor[6 + i] = 0.; // "internal parameter", don't know what for
	}

	// Set up the .csv file for logging
	if (run_counter == 1)
	{
		// Column names
		char* csv_ints_header[] =
		{
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
		};
		char* csv_floats_header[] =
		{
			"TIME",
			"DT",
			"GRAV",
			"RHOW",
			"RHOA",
			"DEPTH",
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
			"timestep_start_time",
			"timestep_end_time",
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
		nr_of_csv_ints = (int)sizeof(csv_ints_header) / sizeof(csv_ints_header[0]);
		nr_of_csv_floats = (int)sizeof(csv_floats_header) / sizeof(csv_floats_header[0]);
		nr_of_csv_doubles = (int)sizeof(csv_doubles_header) / sizeof(csv_doubles_header[0]);
		// Write the log headers
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
	};

	// Log the variables for this timestep
	{
		int csv_ints[] =
		{
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
			iinfo[11]
		};
		float csv_floats[] =
		{
			rinfo[0],
			rinfo[1],
			rinfo[2],
			rinfo[3],
			rinfo[4],
			rinfo[5],
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
			timestep_start_time,
			timestep_end_time,
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
			SAMS_TCP_time, // SAMS_Time
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
			F_ice_global[0],
			F_ice_global[1],
			F_ice_global[2],
			F_ice_global[3],
			F_ice_global[4],
			F_ice_global[5]
		};
		CsvWriter_nextRow(csv_writer);
		char field_buffer[100];
		int field_buffer_length = (int)sizeof(field_buffer) / sizeof(field_buffer[0]);
		for (int i = 0; i < nr_of_csv_ints; i++)
		{
			iResult = sprintf(field_buffer, "%d", csv_ints[i]);
			if (iResult < 0)
			{
				printf("Failure converting integer to .csv field string\n");
			}
			if (iResult >= field_buffer_length)
			{
				printf("Buffer overflow converting integer to .csv field string\n");
			}
			if (CsvWriter_writeField(csv_writer, field_buffer))
			{
				printf("Error writing .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
				*ierr = -23;
				return;
			}
		}
		for (int i = 0; i < nr_of_csv_floats; i++)
		{
			sprintf(field_buffer, "%f", csv_floats[i]);
			if (iResult < 0)
			{
				printf("Failure converting float to .csv field string\n");
			}
			if (iResult >= field_buffer_length)
			{
				printf("Buffer overflow converting float to .csv field string\n");
			}
			if (CsvWriter_writeField(csv_writer, field_buffer))
			{
				printf("Error writing .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
				*ierr = -24;
				return;
			}
		}
		for (int i = 0; i < nr_of_csv_doubles; i++)
		{
			sprintf(field_buffer, "%f", csv_doubles[i]);
			if (iResult < 0)
			{
				printf("Failure converting double to .csv field string\n");
			}
			if (iResult >= field_buffer_length)
			{
				printf("Buffer overflow converting double to .csv field string\n");
			}
			if (CsvWriter_writeField(csv_writer, field_buffer))
			{
				printf("Error writing .csv file: %s\n", CsvWriter_getErrorMessage(csv_writer));
				*ierr = -25;
				return;
			}
		}
	}

	// Last timestep, SAMS won't send anything. Close the connection and logs
	if (substep_nr_overall == nr_of_substeps_overall)
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
		return;
	}

	// Calculations done, prepare for the next timestep
	for (int i = 0; i < 6; i++)
	{
		displacement_SIMA_tm2[i] = displacement_SIMA_tm1[i];
		displacement_SIMA_tm1[i] = displacement_SIMA_t0.double_precision[i];
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

int receive_from_SAMS(SOCKET sams_tcp_socket, double* SAMS_time, double displacement_SAMS[6], double rotation_matrix_SAMS[3][3])
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
	// Rename the Euler axis components and the Euler angle for clarity
	double e1 = SAMS_to_GFEXFO.double_precision[4];
	double e2 = SAMS_to_GFEXFO.double_precision[5];
	double e3 = SAMS_to_GFEXFO.double_precision[6];
	double theta_euler = SAMS_to_GFEXFO.double_precision[7];
	// Check for NaN values
	if (isnan(e1) || isnan(e2) || isnan(e3) || isnan(theta_euler))
	{
		printf("Error at timestep %.1f: axis-angle orientation from SAMS contains bad values\n",*SAMS_time);
		printf("Euler axis: %f %f %f\n", e1, e2, e3);
		printf("Euler angle: %f\n", theta_euler);
		return -1;
	}
	for (int i = 0;i < 3;i++)
		displacement_SAMS[i] = SAMS_to_GFEXFO.double_precision[1 + i]; // Global coordinates
	// Z in SAMS is positive towards the ground, convert here to Z pointing skyward as in SIMA
	e2 = -e2;
	e3 = -e3;
	displacement_SAMS[1] = -displacement_SAMS[1];
	displacement_SAMS[2] = -displacement_SAMS[2];
	
	// Calculations based on data received from SAMS
	rotation_matrix_SAMS[0][0] = (1 - cos(theta_euler)) * e1 * e1 + cos(theta_euler);
	rotation_matrix_SAMS[0][1] = (1 - cos(theta_euler)) * e1 * e2 - e3 * sin(theta_euler);
	rotation_matrix_SAMS[0][2] = (1 - cos(theta_euler)) * e1 * e3 + e2 * sin(theta_euler);
	rotation_matrix_SAMS[1][0] = (1 - cos(theta_euler)) * e2 * e1 + e3 * sin(theta_euler);
	rotation_matrix_SAMS[1][1] = (1 - cos(theta_euler)) * e2 * e2 + cos(theta_euler);
	rotation_matrix_SAMS[1][2] = (1 - cos(theta_euler)) * e2 * e3 - e1 * sin(theta_euler);
	rotation_matrix_SAMS[2][0] = (1 - cos(theta_euler)) * e3 * e1 - e2 * sin(theta_euler);
	rotation_matrix_SAMS[2][1] = (1 - cos(theta_euler)) * e3 * e2 + e1 * sin(theta_euler);
	rotation_matrix_SAMS[2][2] = (1 - cos(theta_euler)) * e3 * e3 + cos(theta_euler);
	bool nan_found = false;
	for (int i = 0;i < 3;i++)
		for (int j = 0;j < 3;j++)
			if (isnan(rotation_matrix_SAMS[i][j]))
				nan_found = true;
	if (nan_found)
	{
		printf("Bad SAMS rotation matrix calculated at timestep %f:\n",*SAMS_time);
		for (int i = 0;i < 3;i++)
		{
			for (int j = 0;j < 3;j++)
				printf("%f ", rotation_matrix_SAMS[i][j]);
			printf("\n");
		}
		printf("Euler axis: %f %f %f\n", e1, e2, e3);
		printf("Euler angle: %f\n", theta_euler);
		return -1;
	}
	
	// Convert the Euler axis-angle representation to Tait-Bryan chained Z-Y-X rotations
	double theta_zyx = asin(-rotation_matrix_SAMS[2][0]);
	double phi_zyx = asin(rotation_matrix_SAMS[2][1] / cos(theta_zyx));
	double psi_zyx = asin(rotation_matrix_SAMS[1][0] / cos(theta_zyx));
	// The calculations were analytically derived from the rotation matrix definition in SIMO manual, app.C
	displacement_SAMS[3] = phi_zyx;
	displacement_SAMS[4] = theta_zyx;
	displacement_SAMS[5] = psi_zyx;
	
	return iResult;
}

int read_from_SAMS(char* SAMS_resultfile_path, int substep_nr_overall, int nr_of_substeps_overall, double* SAMS_txt_time, double F_ice_global[6])
{
	int* ierr = 0;
	int iResult = 0;
	char* line_buffer = NULL; // for getline()
	size_t line_buffer_size = 0; // for getline()
	double SAMS_txt_output[49]; // The last, 49th, column in the .txt file is CalculationTimeRatio [s]
	int header_lines = 7; // The 1st timestep is described by the 8th line in the file
	int lines_to_read = 1;
	static FILE* SAMS_result_file;
	// With a RIFLEX simulation, the first line read will be on sub-timestep 1
	// With a SIMO simulation, the first line read will be on sub-timestep 2
	// Either way, the result file is unopened on the first run
	if (!SAMS_result_file)
	{
		lines_to_read += header_lines;
		printf("Reading SAMS result file: %s\n", SAMS_resultfile_path);
		SAMS_result_file = fopen(SAMS_resultfile_path, "r");
		if (!SAMS_result_file)
		{
			printf("Error opening SAMS result file %s\n",SAMS_resultfile_path);
			*ierr = -2;
			return *ierr;
		}
	}
	int max_wait_ms = 2000;
	static int waited_ms = 0;
	bool starting_over;
	for (int j = 1;j <= max_wait_ms;j *= 5)
	{
		for (int i = 0;i < lines_to_read;i++)
		{
			// TODO handle getline errors everywhere, not just here
			// The fact that a line read failure happens in the middle of the simulation suggests a race condition
			// More investigation is needed, but a possible workaround is to just wait increasing amounts of time
			starting_over = false;
			switch (getline(&line_buffer, &line_buffer_size, SAMS_result_file))
			{
			case 0:
				printf("Empty line i=%d from SAMS result file at sub-timestep %d\n", i, substep_nr_overall);
				break;
			case -1:
				printf("Error %d getting line i=%d from SAMS result file at sub-timestep %d: %s\n", errno, i, substep_nr_overall, strerror(errno));
				*ierr = -18;
				return *ierr;
			case -2:
				starting_over = true;
				if (j > waited_ms) // new waiting record
				{
					waited_ms = j;
					printf("Waited %.3f seconds at sub-timestep %d for line i=%d from %s\n", j / 1000., substep_nr_overall, i, SAMS_resultfile_path);
				}
				Sleep(j);
				lines_to_read = header_lines + substep_nr_overall;
				i = lines_to_read; // Don't attempt to read any more lines
				fseek(SAMS_result_file, 0, SEEK_SET);
				break;
			case -3:
				printf("Line pointer memory allocation failure while getting line i=%d from SAMS result file at sub-timestep %d\n", i, substep_nr_overall);
				*ierr = -31;
				return *ierr;
			case -4:
				printf("Line pointer memory reallocation failure while getting line i=%d from SAMS result file at sub-timestep %d\n", i, substep_nr_overall);
				*ierr = -32;
				return *ierr;
			}
		}
		if (!starting_over) // The switch-case finished successfully, no more waiting or retries needed
			break;
	}
	if (starting_over)
	{
		printf("Timed out (%.3f s) while waiting for last row of %s", max_wait_ms / 1000., SAMS_resultfile_path);
		*ierr = -33;
		return *ierr;
	}
	int current_index = 0;
	int offset = 0;
	// The 29th column is the last ice force component, Ice_Torque_Z_Tot [N.m]
	for (int i = 0;i < 29;i++)
	{
		if (sscanf(line_buffer + current_index, " %lf %n", &SAMS_txt_output[i], &offset) != 1)
		{
			printf("Error parsing column %d from SAMS line:\n%s\n", i, line_buffer);
			*ierr = -19;
			return *ierr;
		}
		current_index += offset;
	}
	// Collect time and ice forces
	*SAMS_txt_time = SAMS_txt_output[0];
	for (int i = 0;i < 3;i++)
	{
		// Linear forces are composed of breaking and rubble forces
		F_ice_global[i] = SAMS_txt_output[20 + i] + SAMS_txt_output[23 + i];
		// Rotational forces are given as total torque in global frame
		F_ice_global[3 + i] = SAMS_txt_output[26 + i];
	}
	return iResult; // 0 if success, negative error code in case of failure
	if (substep_nr_overall == nr_of_substeps_overall)
	{
		fclose(SAMS_result_file); // TODO clean up this file also in case of ungraceful exit
	}
}