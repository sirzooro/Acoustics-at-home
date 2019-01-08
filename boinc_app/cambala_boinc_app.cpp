/*****************************************************************************************
// CAMBALA-based BOINC client application for Acoustics@home -- Copyright(c) 2017
// CAMBALA stands for Coupled Acoustic Modes
// Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory of SB RAS)
// Pavel Petrov (Il'ichev Pacific Oceanological Institute of FEB RAS)
*****************************************************************************************/

#ifdef _WIN32
#include "boinc_win.h"
#else
#include "config.h"
#include <cstdio>
#include <cctype>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <csignal>
#include <unistd.h>
#endif

#include "str_util.h"
#include "util.h"
#include "filesys.h"
#include "boinc_api.h"
#include "mfile.h"

#include <sstream>
#include <string>
#include <fstream>

#include "sequential.h"
#include "utils.h"
#include "point.h"

#if defined(__i386__) || defined (__x86_64__)
#include <cpuid.h>
#endif

#define CHECKPOINT_FILE "chpt"
#define INPUT_FILENAME "in"
#define OUTPUT_FILENAME "out"

__attribute__((noreturn))
void PrintFatalError(const char* str)
{
	// print error to unredirected stderr first
	fprintf(stderr, "Error: %s", str);
	
	// now try to send it back to server
	int retval = boinc_init();
    if (0 == retval)
	{
		fprintf(stderr, "Error: %s", str);
		boinc_finish(1);
    }

	exit(1);
}

void VerifyCpu()
{
#if (defined(__i386__) || defined (__x86_64__)) && defined(__SSE2__)
	unsigned int a, b, c, d;

	if (!__get_cpuid(1, &a, &b, &c, &d))
	{
		PrintFatalError("CPUID instruction is not supported by your CPU!\n");
	}

	if (0 == (d & bit_SSE2))
	{
		PrintFatalError("SSE2 instructions are not supported by your CPU!\n");
	}

#ifdef __SSE4_1__
	if (0 == (c & bit_SSE4_1))
	{
		PrintFatalError("SSE4.1 instructions are not supported by your CPU!\n");
	}
#endif

#ifdef __AVX__
	if (0 == (c & bit_AVX))
	{
		PrintFatalError("AVX instructions are not supported by your CPU!\n");
	}

	// AVX also needs OS support, check for it
	if (0 == (c & bit_OSXSAVE))
	{
		PrintFatalError("OSXSAVE instructions are not supported by your CPU!\n");
	}

	unsigned int eax, edx;
	unsigned int ecx = 0; // _XCR_XFEATURE_ENABLED_MASK
	__asm__ ("xgetbv" : "=a" (eax), "=d" (edx) : "c" (ecx));
	if (0x6 != (eax & 0x6)) // XSTATE_SSE | XSTATE_YMM
	{
		PrintFatalError("AVX instructions are not supported by your OS!\n");
	}
#endif

#ifdef __AVX2__
	if (__get_cpuid_max(0, 0) < 7)
	{
		PrintFatalError("CPUID level 7 is not supported by your CPU!\n");
	}

	unsigned int a2, b2, c2, d2;
	__cpuid_count(7, 0, a2, b2, c2, d2);

	if (0 == (b2 & bit_AVX2))
	{
		PrintFatalError("AVX2 instructions are not supported by your CPU!\n");
	}
#endif

#ifdef __FMA__
	// Some AMD CPUs(s) support FMA but no AVX2. FMA does not provide significat 
	// boost for this app, so it should be enabled on AVX2 CPUs only.
#ifndef __AVX2__
#error AVX2 is not enabled!
#endif
	if (0 == (c & bit_FMA))
	{
		PrintFatalError("FMA instructions are not supported by your CPU!\n");
	}
#endif

#ifdef __AVX512F__
	// AVX512 consists of few subsets. Skylake-AVX512 target supports F, BW, DQ, VL and CD.
	// We need data loaded during AVX and AVX2 checks, make sure we have it
#if !defined(__AVX__) || !defined(__AVX2__)
#error AVX or AVX2 is not enabled!
#endif

	const unsigned int avx512bits = bit_AVX512F | bit_AVX512DQ | bit_AVX512CD | bit_AVX512BW | bit_AVX512VL;

	if (avx512bits != (b2 & avx512bits))
	{
		PrintFatalError("AVX512 instructions are not supported by your CPU!\n");
	}

	// AVX512 also needs OS support, check for it
	if (0xe6 != (eax & 0xe6)) // XSTATE_SSE | XSTATE_YMM | XSTATE_OPMASK | XSTATE_ZMM | XSTATE_HI_ZMM
	{
		PrintFatalError("AVX512 instructions are not supported by your OS!\n");
	}
#endif

#endif
}

using namespace std;
char buf[256];

int do_work(search_space_point &cur_record_point);
int do_checkpoint( const long long &total_points, 
	               const long long &processed_points, 
	               const search_space_point &cur_record_point );

int main(int argc, char **argv)
{
	VerifyCpu();

	int retval = boinc_init();
    if ( retval ) {
        fprintf(stderr, "%s APP: boinc_init() returned %d\n",
            boinc_msg_prefix(buf, sizeof(buf)), retval
        );
        exit( retval );
    }

	search_space_point cur_record_point;
	cur_record_point.residual = START_HUGE_VALUE;

	retval = do_work(cur_record_point);
	if ( retval ) {
		fprintf( stderr, "APP: do_work() returned %d \n", retval);
		exit(retval);
	}

	// resolve, open and write answer to output file
	string output_file_name = "";
    retval = boinc_resolve_filename_s( OUTPUT_FILENAME, output_file_name);
	if (retval) {
		fprintf(stderr, "%s APP: boinc_resolve_filename_s() returned %d \n",
			boinc_msg_prefix(buf, sizeof(buf)), retval
		);
		exit(retval);
	}
	if (output_file_name == "") {
		fprintf(stderr, "APP: output_file_name is empty \n");
		exit(-1);
	}

	ofstream output_file(output_file_name.c_str());
    if ( !output_file.is_open() ) {
        fprintf(stderr, "APP: output file opening failed \n");
		exit(-1);
    }
	
	if (!isCorrectCalculatedPoint(cur_record_point)) {
		fprintf(stderr, "APP: incorrect final point to write \n");
		return -1;
	}
	fromPointToFile(cur_record_point, output_file);
	
	output_file.close();
	
    boinc_finish(0);
}

int do_work (search_space_point &cur_record_point)
{
	string input_file_name = "";
	string str = "";

	// open the input file (resolve logical name first)
	int retval = boinc_resolve_filename_s(INPUT_FILENAME, input_file_name);
	if (retval) {
		fprintf(stderr, "%s APP: boinc_resolve_filename_s() returned %d \n",
			boinc_msg_prefix(buf, sizeof(buf)), retval
		);
		exit(retval);
	}
	if (input_file_name == "") {
		fprintf(stderr, "APP: input_file_name is empty \n");
		exit(-1);
	}

	CAMBALA_sequential cambala_seq;
	retval = cambala_seq.readScenario(input_file_name);
	if (retval) {
		fprintf(stderr, "APP: readScenario() returned %d \n", retval);
		exit(retval);
	}

	retval = cambala_seq.readInputDataFromFiles();
	if (retval) {
		fprintf(stderr, "APP: readInputDataFromFiles() returned %d \n", retval);
		exit(retval);
	}
	
	vector<vector<double>> depths_vec = cambala_seq.createDepthsArray();
	if (depths_vec.size() > 1) {
		fprintf(stderr, "APP: depths_vec.size() %zd \n", depths_vec.size());
		exit(-1);
	}
	vector<double> depths = depths_vec[0];
	depths_vec.clear();
	cout << "depths :";
	for (unsigned i = 0; i < depths.size(); i++)
		cout << " " << depths[i];
	cout << endl;
	retval = cambala_seq.init(depths);
	if (retval) {
		fprintf(stderr, "APP: init() failed %d\n", retval);
		exit(retval);
	}
	
	vector<search_space_point> points_vec = cambala_seq.getSearchSpacePointsVec(depths);
	long long total_points = points_vec.size();
	
	if (total_points <= 0) {
		fprintf(stderr, "APP: do_checkpoint() total_points <= 0 \n");
		return -1;
	}

	// read data from the checkpoint file if such exists
	string chpt_file_name = "";
	long long processed_points = 0;
	retval = boinc_resolve_filename_s(CHECKPOINT_FILE, chpt_file_name);
	if (retval) {
		fprintf(stderr, "%s APP: in do_work() boinc_resolve_filename_s() returned %d \n",
			boinc_msg_prefix(buf, sizeof(buf)), retval
		);
		exit(retval);
	}
	if (chpt_file_name == "") {
		fprintf(stderr, "APP: in do_work() chpt_file_name is empty \n");
		exit(-1);
	}
	ifstream chpt_file;
	chpt_file.open(chpt_file_name.c_str(), ios_base::in);
	if (chpt_file.is_open()) {
		str = "";
		getline(chpt_file, str);
		if (str == "") {
			fprintf(stderr, "APP: in do_work() couldn't read processed_points from the chpt file \n");
			return -1;
		}
		istringstream(str) >> processed_points;
		if ( (processed_points <= 0) || (processed_points > total_points) ) {
			fprintf(stderr, "APP: in do_work() incorrect processed_points \n");
			return -1;
		}
		str = "";
		getline(chpt_file, str);
		if (str == "") {
			fprintf(stderr, "APP: in do_work() couldn't read point from the chpt file \n");
			return -1;
		}
		unsigned long long cws_count = cambala_seq.cw1_arr.size();
		if (!cws_count) {
			fprintf(stderr, "APP: in do_work() cws_count == 0 \n");
			return -1;
		}
		cur_record_point = fromStrToPoint(str, cws_count);
		if (!isCorrectCalculatedPoint(cur_record_point)) {
			fprintf(stderr, "APP: in do_work() incorrect final point to write \n");
			return -1;
		}
		cambala_seq.updateRecordPoint(cur_record_point);
		chpt_file.close();
		cout << "point from chpt file" << endl;
	}
	if (processed_points == total_points) // exit if all points were already processed
		return 0;
	
	double dval = 0;
	long long checkpoint_every_point = 1;
	if ((cambala_seq.object_function_type == "uniform") || (cambala_seq.object_function_type == "weighted"))
		checkpoint_every_point = 10;
	
	for ( long long i = processed_points; i < total_points; i++) {
		dval = cambala_seq.fillDataComputeResidual(points_vec[i]);
		if (dval <= 0) {
			fprintf(stderr, "APP: in do_work() residual <= 0 \n");
			return -1;
		}
		if (dval < cur_record_point.residual) {
			cur_record_point = points_vec[i];
			if (!isCorrectCalculatedPoint(cur_record_point)) {
				fprintf(stderr, "APP: in do_work() incorrect point \n");
				return -1;
			}
		}
		
		// checkpoint current results
		//if ( ( boinc_is_standalone() ) || ( boinc_time_to_checkpoint() ) ) {
		if ((i+1) % checkpoint_every_point == 0) {
			retval = do_checkpoint(total_points, i + 1, cur_record_point);
			if (retval) {
				fprintf(stderr, "APP: checkpoint failed %d \n", retval);
				exit(retval);
			}
			boinc_checkpoint_completed();
		}
		cout << "processed " << i+1 << " out of " << total_points << endl;
        //}
	}
	
	return 0;
}

int do_checkpoint( const long long &total_points, 
	               const long long &processed_points, 
	               const search_space_point &cur_record_point )
{
	//check input data 
	if (processed_points <= 0) {
		fprintf(stderr, "APP: do_checkpoint() processed_points <= 0 \n");
		return -1;
	}

	if (total_points <= 0) {
		fprintf(stderr, "APP: do_checkpoint() total_points <= 0 \n");
		return -1;
	}

	if (!isCorrectCalculatedPoint(cur_record_point)) {
		fprintf(stderr, "APP: in do_checkpoint() an incorrect point to write \n");
		return -1;
	}
	
	ofstream temp_ofile("temp");
	if (!temp_ofile.is_open()) {
		fprintf(stderr, "APP: in do_checkpoint() temp ofile wasn't opened \n");
		return -1;
	}
	
	temp_ofile << processed_points << endl;
	fromPointToFile(cur_record_point, temp_ofile);
    temp_ofile.close();
	
	string resolved_name = "";
    boinc_resolve_filename_s(CHECKPOINT_FILE, resolved_name);
	if (resolved_name == "") {
		fprintf(stderr, "APP: resolved_name is empty \n");
		return -1;
	}
	
	int retval = 0;
    retval = boinc_rename( "temp", resolved_name.c_str() );
	if ( retval ) {
		fprintf(stderr, "APP: do_checkpoint() boinc_rename() returned %d \n", retval);
		return retval;
	}

	boinc_fraction_done( (double)processed_points / (double)total_points );

    return 0;
}

#ifdef _WIN32
int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR Args, int WinMode) {
    LPSTR command_line;
    char* argv[100];
    int argc;

    command_line = GetCommandLine();
    argc = parse_command_line( command_line, argv );
    return main(argc, argv);
}
#endif