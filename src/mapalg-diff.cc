/** mapalg-diff.cc: map algebra - diff two rasters
 * Author: Yan Y. Liu <yanliu@illinois.edu>
 * Date: 06/18/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
#include "data.h"
#include "gdal.h"
#include "cpl_conv.h"
#include "cpl_string.h"

#include "mpi.h"

// assumption: float value; same dimensions
int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	int rank, np;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	double t0, t1, t2, t3, t4, t5; // timing info

	t0 = get_timemark();
	// step 1: read input rasters
	char * fn1 = argv[1]; // raster 1
	char * fn2 = argv[2]; // raster 2
	char * ofn = argv[3]; // output: diff raster
	int i, j;

	// data only proc 0 knows and handles	
	double georef[6]; // georef data structure for a raster
	char prj[2048]; // store projection wkt
	double nodata;
	int x=0, y=0; // size of raster on x and y dim
	float *raster1 = NULL, *raster2 = NULL;
	int *offsetx = NULL, *offsety = NULL, *sizex = NULL, *sizey = NULL;
	// data each proc needs to know
	int maxx=0, maxy=0, boffsetx = 0, boffsety = 0, bsizex = 0, bsizey = 0;
	float *b1, *b2; // two blocks
	float *block; // tmp pointer in accessing global raster to a block

    if (rank == 0) { // root proc reads dataset
	GDALDatasetH rin1, rin2;
	// get input raster info
	rin1 = raster_open(fn1, georef, prj, &nodata, &x, &y);
	raster_info(rin1, georef, prj, nodata, x, y);
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
	raster1 = (float *)malloc(sizeof(float) * maxx * maxy * np);
	memset(raster1, 0, sizeof(float) * maxx * maxy * np);
	// read blocks from input raster
	for (i=0; i<np; i++) {
		block = raster1 + i * (maxx * maxy);
		raster_read(rin1, block, offsetx[i], offsety[i], sizex[i], sizey[i]);
	}
	// close raster
	raster_close(rin1);
	// read 2nd raster
	rin2 = raster_open(fn2, georef, prj, &nodata, &x, &y);
        raster_info(rin2, georef, prj, nodata, x, y); // they should be the same
	raster2 = (float *)malloc(sizeof(float) * maxx * maxy * np);
	memset(raster2, 0, sizeof(float) * maxx * maxy * np);
	// read blocks from input raster
	for (i=0; i<np; i++) {
		block = raster2 + i * (maxx * maxy);
		raster_read(rin2, block, offsetx[i], offsety[i], sizex[i], sizey[i]);
	}
	// close raster
	raster_close(rin2);
    }

	t1 = get_timemark();
	// step 2: transfer data blocks to procs
	// scatter data blocks to procs
	MPI_Bcast(&maxx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&maxy, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nodata, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(offsetx, 1, MPI_INT, &boffsetx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(offsety, 1, MPI_INT, &boffsety, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(sizex, 1, MPI_INT, &bsizex, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(sizey, 1, MPI_INT, &bsizey, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG
	fprintf(stderr, "rank %d: max[x,y]=%d,%d nodata=%.5lf offset=%d,%d size=%d,%d\n", rank, maxx, maxy, nodata, boffsetx, boffsety, bsizex, bsizey);
#endif
	// allocate memory for blocks
	b1 = (float *) malloc(sizeof(float) * maxx * maxy);
	b2 = (float *) malloc(sizeof(float) * maxx * maxy);
	MPI_Scatter(raster1, maxx * maxy, MPI_FLOAT, b1, maxx * maxy, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Scatter(raster2, maxx * maxy, MPI_FLOAT, b2, maxx * maxy, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	t2 = get_timemark();
	// step 3: map algebra operation - diff
	int index;
	for (i=0; i<bsizey; i++) {
		for (j=0; j<bsizex; j++) {
			index = i * bsizex + j;
			if (fabs((double)(b1[index]) - nodata) < 0.001)
				b1[index] = 0;
			else 
				b1[index] = b1[index] - b2[index];
		}
	}
	
	t3 = get_timemark();
	// step 4: transfer results for writing
	MPI_Gather(b1, maxx * maxy, MPI_FLOAT, raster1, maxx * maxy, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	t4 = get_timemark();
	// step 5: write output
    if (rank == 0) {
	GDALDatasetH rout;
	rout = raster_create(ofn, x, y, georef, prj, 0, 0);
	for (i=0; i<np; i++) {
		block = raster1 + i * (maxx * maxy);
		raster_write(rout, block, offsetx[i], offsety[i], sizex[i], sizey[i]);
	}
	raster_close(rout);
    }

	t5 = get_timemark();
    if (rank == 0) {
	jobstat.Tread = t1 - t0;
	jobstat.Tcommdata = t2 - t1;
	jobstat.Tcompute = t3 - t2;
	jobstat.Tcommresult = t4 - t3;
	jobstat.Twrite = t5 - t4;
	jobstat.Ttotal = t5 - t0;
	print_jobstat();
    }
	// step 6: clean up
	free(b1);
	free(b2);
    if (rank == 0) {
	free(offsetx);
	free(offsety);
	free(sizex);
	free(sizey);
	free(raster1);
	free(raster2);
    }
	MPI_Finalize();
}
