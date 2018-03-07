/** mapalg-dist.cc: map algebra - focal distance from a point
 * multi read single write. no MPI IO though
 * Author: Yan Y. Liu <yanliu@illinois.edu>
 * 06/20/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
#include "data.h"
#include "gdal.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include "ogr_spatialref.h"

#include "mpi.h"

#include "sptw.h"

// assumption: float value on DEM; calc dist using 3 dims
// command: mapalg-dist dem px py pz dist_raster
// e.g.: mpirun -np 4 mapalg-dist ../../data/allbig.tif -1492833 2732555 1170 ../../output/dist.tif
int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	int rank, np;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	double t0, t1, t2, t3, t4, t5; // timing info

	t0 = get_timemark();
	// step 1: read input rasters
	char * fn = argv[1]; // DEM
	float xp = (float)atof(argv[2]); // x  
	float yp = (float)atof(argv[3]); // y  
	float zp = (float)atof(argv[4]); // z  
	char * ofn = argv[5]; // output: dist raster
	int i, j;

	// data: every proc knows and handles	
	double georef[6]; // georef data structure for a raster
	char prj[2048]; // store projection wkt
	double nodata;
	int x=0, y=0; // size of raster on x and y dim
	float *raster = NULL;
	int *offsetx = NULL, *offsety = NULL, *sizex = NULL, *sizey = NULL;
	// data each proc needs to know
	int maxx=0, maxy=0, boffsetx = 0, boffsety = 0, bsizex = 0, bsizey = 0;
	float *b; // two blocks
	float *block; // tmp pointer in accessing global raster to a block

	// every proc reads its own block
	GDALDatasetH rin;
	// get input raster info
	rin = raster_open(fn, georef, prj, &nodata, &x, &y);
	raster_info(rin, georef, prj, nodata, x, y);
	// determine block sizes
	get_block(rank, np, x, y, &boffsetx, &boffsety, &bsizex, &bsizey);
	// we still need to find max sizex and sizey for output, do it MPI way
	MPI_Allreduce(&bsizex, &maxx, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&bsizey, &maxy, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	// allocate data block memory
	block = (float *)malloc(sizeof(float) * maxx * maxy);
	memset(block, 0, sizeof(float) * maxx * maxy);
	// read blocks from input raster
	raster_read(rin, block, boffsetx, boffsety, bsizex, bsizey);
	// close raster
	raster_close(rin);

	t1 = get_timemark();
	// step 2: transfer data blocks to procs
	// no need to communicate metadata since everyone has it
#ifdef DEBUG
	fprintf(stderr, "rank %d: max[x,y]=%d,%d nodata=%.5lf offset=%d,%d size=%d,%d\n", rank, maxx, maxy, nodata, boffsetx, boffsety, bsizex, bsizey);
#endif

	t2 = get_timemark();
	// step 3: map algebra operation - focal distance
	b = block;
	int index;
	float xc, yc, zc;
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
	
	t3 = get_timemark();
	// step 4: transfer results for writing
	// since each proc writes, no need to transfer output

	t4 = get_timemark();
	// step 5: write output: every proc writes
    if (rank == 0) {
	// create output raster
	if (sptw::create_tiled_raster(ofn, x, y, 1, GDT_Float32, georef, prj, 4096) != sptw::SP_None) {
		fprintf(stderr, "ERROR: create output raster %s\n", ofn);
		exit(1);
	}
#ifdef DEBUG
	fprintf(stderr, "raster file %s created\n", ofn);
#endif
    }
	MPI_Barrier(MPI_COMM_WORLD);
	// init output raster
	sptw::PTIFF* rtmp = sptw::open_raster(ofn);
    if (rank == 0) {
	if (sptw::populate_tile_offsets(rtmp, 4096) != sptw::SP_None) {
		fprintf(stderr, "ERROR: create init offset for geotiff %s\n", ofn);
		exit(1);
	}
#ifdef DEBUG
	fprintf(stderr, "raster file %s offset inited\n", ofn);
#endif
    }
	sptw::close_raster(rtmp); // TODO for sptw: move this to create func
	MPI_Barrier(MPI_COMM_WORLD);
	
	// open output raster for writing by each proc
	sptw::PTIFF* rout = sptw::open_raster(ofn); 
	//sptw::write_area(rout, block, (long int)boffsetx, (long int)boffsety, (long int)(boffsetx + sizex - 1), (long int)(boffsety + sizey -1));
	sptw::write_area(rout, block, (long int)boffsetx, (long int)boffsety, (long int)(boffsetx + 4096 - 1), (long int)(boffsety + 4096 - 1));
	MPI_Barrier(MPI_COMM_WORLD);
	sptw::close_raster(rout);
	MPI_Barrier(MPI_COMM_WORLD);

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
	free(b);
    if (rank == 0) {
	free(offsetx);
	free(offsety);
	free(sizex);
	free(sizey);
	free(raster);
    }
	MPI_Finalize();
}
