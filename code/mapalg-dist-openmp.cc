/** mapalg-dist.cc: map algebra - focal distance from a point
 * using openmpi parallel programming model.
 * Author: Yan Y. Liu <yanliu@illinois.edu>
 * Date: 06/18/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "util.h"
#include "data.h"
#include "gdal.h"
#include "cpl_conv.h"
#include "cpl_string.h"

// assumption: float value
// command: mapalg-dist dem px py pz dist_raster num_threads
// e.g.: export OMP_NUM_THREADS=32; mapalg-dist-openmp ../../data/allbig.tif -1492833 2732555 1170 ../../output/dist.tif 4
int main(int argc, char** argv) {
	int np = atoi(argv[6]); 
	omp_set_num_threads (np); // set by env var OMP_NUM_THREADS
	double t0, t1, t2, t3, t4, t5; // timing info

	t0 = get_timemark();
	// step 1: read input rasters
	char * fn = argv[1]; // DEM
	float xp = (float)atof(argv[2]); // x  
	float yp = (float)atof(argv[3]); // y  
	float zp = (float)atof(argv[4]); // z  
	char * ofn = argv[5]; // output: dist raster
	int i, j;

	// data only proc 0 knows and handles	
	double georef[6]; // georef data structure for a raster
	char prj[2048]; // store projection wkt
	double nodata;
	int x=0, y=0; // size of raster on x and y dim
	float *raster = NULL;
	int *offsetx = NULL, *offsety = NULL, *sizex = NULL, *sizey = NULL;
	// data each thread needs to know
	int maxx=0, maxy=0, boffsetx = 0, boffsety = 0, bsizex = 0, bsizey = 0;
	float *b; // two blocks
	float *block; // tmp pointer in accessing global raster to a block

	GDALDatasetH rin;
	// get input raster info
	rin = raster_open(fn, georef, prj, &nodata, &x, &y);
	raster_info(rin, georef, prj, nodata, x, y);
	// determine block sizes
	offsetx = (int *) malloc(sizeof(int) * np);
	memset(offsetx, 0, sizeof(int) * np);
	offsety = (int *) malloc(sizeof(int) * np);
	memset(offsety, 0, sizeof(int) * np);
	sizex = (int *) malloc(sizeof(int) * np);
	memset(sizex, 0, sizeof(int) * np);
	sizey = (int *) malloc(sizeof(int) * np);
	memset(sizey, 0, sizeof(int) * np);
	for (i=0; i<np; i++) {
		get_block(i, np, x, y, &offsetx[i], &offsety[i], &sizex[i], &sizey[i]);
	}
	// find max sizex and sizey
	for (i=0; i<np; i++) {
		if (maxx < sizex[i]) maxx = sizex[i];
		if (maxy < sizey[i]) maxy = sizey[i];
	}
	// allocate data block memory
	raster = (float *)malloc(sizeof(float) * maxx * maxy * np);
	memset(raster, 0, sizeof(float) * maxx * maxy * np);
	// read blocks from input raster
	for (i=0; i<np; i++) {
		block = raster + i * (maxx * maxy);
		raster_read(rin, block, offsetx[i], offsety[i], sizex[i], sizey[i]);
	}
	// close raster
	raster_close(rin);

	t1 = get_timemark();
	// step 2: transfer data blocks to procs

	t2 = get_timemark();
	// step 3: map algebra operation - focal distance
	int index;
	float xc, yc, zc;
	int p;
	#pragma omp parallel
	fprintf(stdout, "OpenMP mode: using %d threads\n", omp_get_num_threads()); 
	#pragma omp parallel for private(p, i, j, index, xc, yc, zc, boffsetx, boffsety, bsizex, bsizey, b)
	for (p=0; p<np; p++) {
		boffsetx = offsetx[p];
		boffsety = offsety[p];
		bsizex = sizex[p];
		bsizey = sizey[p];
		b = raster + p * (maxx * maxy);
		for (i=0; i<bsizey; i++) {
			for (j=0; j<bsizex; j++) {
				index = i * bsizex + j;
				if (b[index] - nodata < 0.001) continue;
				// get coord of cell. TODO: use centroid
				xc = georef[0] + (boffsetx + j) * georef[1];
				yc = georef[3] + (boffsety + i) * georef[5];
				zc = b[index];
				b[index] = get_dist(xc, yc, zc, xp, yp, zp);
			}
		}
	}
	t3 = get_timemark();
	// step 4: transfer results for writing

	t4 = get_timemark();
	// step 5: write output
	GDALDatasetH rout;
	rout = raster_create(ofn, x, y, georef, prj, nodata);
	for (i=0; i<np; i++) {
		block = raster + i * (maxx * maxy);
		raster_write(rout, block, offsetx[i], offsety[i], sizex[i], sizey[i]);
	}
	raster_close(rout);

	t5 = get_timemark();
	jobstat.Tread = t1 - t0;
	jobstat.Tcommdata = t2 - t1;
	jobstat.Tcompute = t3 - t2;
	jobstat.Tcommresult = t4 - t3;
	jobstat.Twrite = t5 - t4;
	jobstat.Ttotal = t5 - t0;
	print_jobstat();

	// step 6: clean up
	free(offsetx);
	free(offsety);
	free(sizex);
	free(sizey);
	free(raster);
}
