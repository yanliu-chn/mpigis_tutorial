#ifndef UTIL_H
#define UTIL_H
/** util.h: common utility functions needed by other functions
 * Author: Yan Y. Liu <yanliu@illinois.edu>
 * Date: 06/18/2014
 */

typedef struct {
	double Tread;
	double Tcommdata;
	double Tcompute;
	double Tcommresult;
	double Twrite;
	double Ttotal;
} Jobstat;
extern Jobstat jobstat;

float get_dist(float x1, float y1, float z1, float x2, float y2, float z2);
double get_timemark();
void print_jobstat() ;
int get_best_dim(int np, int *rowSize, int *colSize);
int get_block(int rank, int np, int x, int y, int *offsetx, int *offsety, int * sizex, int *sizey);
int get_block2(int rank, int np, int x, int y, int *offsetx, int *offsety, int * sizex, int *sizey);

#endif
