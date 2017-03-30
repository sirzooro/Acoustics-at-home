#include "bisect_cpu.h"
#include "assert.h"
#include "stdlib.h"
#include <iostream>
#include <math.h>


static const float precision = 0.000001;
#define MIN_ABS_INTERVAL (float)5.0e-37

int computeSmallerEigenvals (const float x, const float* md, const float* sd, const int mat_size)
{
	//assert (sd[0]==0.0f);
	// perform (optimized) Gaussian elimination to determine the number
	// of eigenvalues that are smaller than x
	const float* sd_shifted = sd-1;
	float  delta = 1.0f;
	int count = 0;
	delta = md[0] - x;
	count += (delta < 0);
	for (int i = 1; i < mat_size; ++i)
	{
		delta = md[i] - x - (sd_shifted[i] * sd_shifted[i]) / delta;
		count += (delta < 0);
	}
	assert (count>=0);
	assert (count<=mat_size);
	return count;

}

float midpoint(const float& left, const float& right)
{
	// try to avoid overflow
	float mid;
	if (sign_f (left) == sign_f (right))
		mid = left + (right - left) * 0.5f;
	else
		mid = (left + right) * 0.5f;

	return mid;
}

bool intervalConverged (const Interval& c)
{
	float t0 = fabs (c.rl - c.ll);
	float t1 = max (fabs(c.ll), fabs(c.rl)) * precision;
	return (t0 <= max (MIN_ABS_INTERVAL, t1));

}


void PushOrStoreInterval (const Interval curint, const Interval_e curint_e,
		Interval* intervals, Interval_e* intervals_e, int& intervals_c,
		float* eigenvalues, int& eigenvalues_c)
{
	if (intervalConverged (curint))
		for (int i=0; i<(curint_e.rc-curint_e.lc); ++i)
			StoreEv (midpoint (curint), eigenvalues, eigenvalues_c);
	else
		PushInterval (curint, curint_e, intervals, intervals_e, intervals_c);
	//std::cout <<  "interval  " <<  curint.ll << " " << curint.rl << " " << intervalConverged(curint) << "\n" << std::flush;
}

int bisectCpu (const float* md, const float* sd, const int mat_size,
		const float ll, const float rl,
		float* eigenvalues)
{
	assert(ll <= rl);
	int eigenvalues_c = 0; // Number of eigenvalues found

	// We can't put these into a single structure, since CUDA will
	// drop it from registers into local memory
	Interval intervals[mat_size]; // Intervals borders
	Interval_e intervals_e[mat_size]; // Eigenvalue counts for intervals borders
	int intervals_c = 0; // Intervals stack height

	Interval_e toplevel_e = Interval_e {
				computeSmallerEigenvals (ll, md, sd, mat_size),
				computeSmallerEigenvals (rl, md, sd, mat_size)};
	PushOrStoreInterval (Interval {ll, rl}, toplevel_e,  intervals, intervals_e, intervals_c, eigenvalues, eigenvalues_c);

	while (intervals_c > 0)
	{
		// Pop intervals stack
		--intervals_c;
		Interval_e curint_e = intervals_e[intervals_c];
		Interval curint = intervals[intervals_c];
		int lc = curint_e.lc;
		int rc = curint_e.rc;
		//std::cout << mat_size << " " << intervals_c << " Lim: " << curint.ll << " " << curint.rl << "lrc " << lc << " " << rc << "\n" << std::flush;

		float m = midpoint (curint);
		int mc = computeSmallerEigenvals (m, md, sd, mat_size);
		if (mc-lc > 0)
			PushOrStoreInterval (Interval {curint.ll, m}, Interval_e {lc, mc}, intervals, intervals_e, intervals_c, eigenvalues, eigenvalues_c);
		if (rc-mc > 0)
			PushOrStoreInterval (Interval {m, curint.rl}, Interval_e {mc, rc}, intervals, intervals_e, intervals_c, eigenvalues, eigenvalues_c);
	}
	return eigenvalues_c;
}


