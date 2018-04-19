#ifndef UTIL_CC
#define UTIL_CC
/** util.cc: common utility functions needed by other functions
 * Author: Yan Y. Liu <yanliu@illinois.edu>
 * Date: 06/18/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "util.h"

Jobstat jobstat; // stat 

// distance calcuation: considers elevation. must have the same unit for all
// coordinates
float get_dist(float x1, float y1, float z1, float x2, float y2, float z2)
{
	return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}

// get current system time
double get_timemark()
{
	struct timeval tsec;
	struct timezone tzone;
	gettimeofday(&tsec, &tzone);
	return (double)(tsec.tv_sec + tsec.tv_usec/1000000.0);
}

void print_jobstat() 
{
	fprintf(stdout, "====================\n");
	fprintf(stdout, "||    Job Stat    ||\n");
	fprintf(stdout, "====================\n");
	fprintf(stdout, "Reading time: %.5lf seconds\n", jobstat.Tread);
	fprintf(stdout, "Data distribution time: %.5lf seconds\n", jobstat.Tcommdata);
	fprintf(stdout, "Running time: %.5lf seconds\n", jobstat.Tcompute);
	fprintf(stdout, "Result collection time: %.5lf seconds\n", jobstat.Tcommresult);
	fprintf(stdout, "Writing time: %.5lf seconds\n", jobstat.Twrite);
	fprintf(stdout, "TOTAL TIME: %.5lf seconds\n", jobstat.Ttotal);
}

// find the rectangle size closest to the square root of np
int get_best_dim(int np, int *rowSize, int *colSize)
{
	int i = 0; // loop counter
	int p = (int) sqrt(np);
	int q = np / p;
	while (p > 1 &&(p * q != np || p > q)) {
		p --;
		q = np / p;
		i++;
	}
	*rowSize = p;
	*colSize = q;
	return i;
}

// calc the data block to be read given the rank of the process
// rank, np: rank and number of procs
// x, y: size of raster 
// offsetx, offsety: output. 
// sizex, sizey: output. useful for procs positioned at the end of each dim
int get_block(int rank, int np, int x, int y, int *offsetx, int *offsety, int * sizex, int *sizey)
{
	int dimx, dimy, blockx, blocky;
	int cellsx, cellsy, rcellsx, rcellsy; // num cells on 2 dims and remainders
	get_best_dim(np, &dimy, &dimx); // number of blocks on each dimension
	blockx = rank % dimx; // block id on x dim
	blocky = rank / dimx; // block id on y dim
#ifdef DEBUG
	fprintf(stderr, "np %d gets %d x %d blocks; rank %d blocky=%d blockx=%d\n", np, dimy, dimx, rank, blocky, blockx);
#endif
	cellsx = x / dimx;
	cellsy = y / dimy;
	rcellsx = x % dimx;
	rcellsy = y % dimy;
	*offsetx = blockx * cellsx;
	*offsety = blocky * cellsy;
	*sizex = cellsx;
	*sizey = cellsy;
	if (blockx == dimx - 1) { // end block on x dim
		*sizex = *sizex + rcellsx;
	}
	if (blocky == dimy - 1) { // end block on y dim
		*sizey = *sizey + rcellsy;
	}
	return 1;
}


int get_block2(int rank, int np, int x, int y, int *offsetx, int *offsety, int * sizex, int *sizey)
{
	int dimx, dimy, blockx, blocky;
	int cellsx, cellsy, rcellsx, rcellsy; // num cells on 2 dims and remainders
	get_best_dim(np, &dimy, &dimx); // number of blocks on each dimension
	blockx = rank % dimx; // block id on x dim
	blocky = rank / dimx; // block id on y dim
#ifdef DEBUG
	fprintf(stderr, "np %d gets %d x %d blocks; rank %d blocky=%d blockx=%d\n", np, dimy, dimx, rank, blocky, blockx);
#endif
	cellsx = x / dimx;
	cellsy = y / dimy;
	rcellsx = x % dimx;
	rcellsy = y % dimy;
	*offsetx = blockx * cellsx;
	*offsety = blocky * cellsy;
	*sizex = cellsx;
	*sizey = cellsy;
	if (blockx < rcellsx) { *offsetx = *offsetx + rank; *sizex = *sizex + 1; }
	else { *offsetx = *offsetx + rcellsx; }
	if (blocky < rcellsy) { *offsety = *offsety + rank; *sizey = *sizey + 1; }
	else { *offsety = *offsety + rcellsy; }
	return 1;
}
#endif
