/*****************************************************************************************
// CAMBALA: Coupled Acoustic Modes -- Copyright(c) 2015-2017
// Pavel Petrov (Il'ichev Pacific Oceanological Institute of FEB RAS), 
// Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory of SB RAS)
*****************************************************************************************/

#ifdef _MPI
#include <mpi.h>
#endif

#include "sequential.h"
#include "parallel.h"

#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cmath>

int main(int argc, char **argv)
{
	unsigned ncpl = 0; // search mesh within each water layer

	std::string scenarioFileName = "";
	CAMBALA_sequential CAMBALA_seq;
	int verbosity = 0;

#ifdef _DEBUG
	argc = 3;
	scenarioFileName = "..//scenarios//test_fixed_depths.txt";
	//argv[1] = "./boinc_CAMBALA_app/in";
	//argv[1] = "39_hydro_r_uniform260.txt";
	//argv[1] = "true_scenario_2.txt";
	//argv[2] = "1";
	verbosity = 2;
#else
	if (argc > 1)
		scenarioFileName = argv[1];
	else {
		cout << "Usage : scenarioFileName [verbosity] [-test]" << endl;
		exit(1);
	}
	if (argc >= 3)
		verbosity = atoi(argv[2]);
#endif
	
#ifndef _MPI
	// sequential mode%
	cout << "verbosity " << verbosity << endl;

	// fix start time
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	std::chrono::high_resolution_clock::time_point t2;
	std::chrono::duration<double> time_span;

	CAMBALA_seq.verbosity = verbosity;
	// read scenario, modal_delays, mode_numbers and freqs, then determine the search space
	CAMBALA_seq.readScenario(scenarioFileName);
	CAMBALA_seq.readInputDataFromFiles();
	
	CAMBALA_seq.solve();
	CAMBALA_seq.reportFinalResult();

	t2 = std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	cout << "main() total time " << time_span.count() << " s" << endl;

#else
	int rank = 0;
	int corecount = 1;

	// parallel mode
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &corecount);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (!rank)
		cout << "verbosity " << verbosity << endl;

	CAMBALA_parallel CAMBALA_par;
	CAMBALA_par.rank = rank;
	CAMBALA_par.corecount = corecount;
	CAMBALA_par.verbosity = verbosity;

	CAMBALA_par.readScenario(scenarioFileName);
	CAMBALA_par.readInputDataFromFiles();
	
	double cur_time = MPI_Wtime();
	CAMBALA_par.MPI_main();
	cout << "MPI_main() total time " << MPI_Wtime() - cur_time << " s" << endl;
#endif
	
	return 0;
}