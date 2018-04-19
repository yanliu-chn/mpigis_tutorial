/** rastercp.cc: Illustration of GDAL C API-based raster copy
 * Author: Yan Y. Liu <yanliu@illinois.edu>
 * Date: 06/18/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "util.h"
#include "data.h"
#include "gdal.h"
#include "cpl_conv.h"
#include "cpl_string.h"

// create test raster file 
// g++ -DDEBUG -I. -I../code/ -I$GDAL_HOME/include -o gentestdata gentestdata.cc ../code/util.cc ../code/data.cc -L$GDAL_HOME/lib -lgdal -lm
int main(int argc, char **argv) {
	char * ofn = argv[1];
	int i, j;
	
	double georef[6] = {-80.0, 0.001, 0, 40.0, 0, -0.01}; // georef data structure for a raster
	char prj[2048]; // store projection wkt
	double nodata = 0;
	int x=10000, y=10000; // size of raster on x and y dim
	OGRSpatialReference srs;
	srs.SetFromUserInput("EPSG:4326");
	srs.exportToPrettyWkt(&prj);

	// allocate data block memory
	float *raster = (float *)malloc(sizeof(float) * x * y);
	memset(raster, 0, sizeof(float) * x * y );
	// read blocks from input raster
	int I = 10;
	for (i=0; i<y; i++) {
		for (j=0; j<x; j++) {
			if (j % I == i % I)
				raster[i * x + j] = 0;
			else
				raster[i * x + j] = i * x + j;
		}
	}
	// create output raster
	GDALDatasetH rout;
	rout = raster_create(ofn, x, y, georef, prj, nodata);
	raster_write(rout, raster, 0, 0, x, y);
	raster_close(rout);
	
	// free memory
	free(raster);
}
