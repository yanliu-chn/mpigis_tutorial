/** sendrecv.cc: basic send-receive program to illustrate point2point comm
 * and deadlock situations
 * Only works if np > 0 and np is an even number
 * Author: Yan Y. Liu <yanliu@illinois.edu>
 * Date: 06/18/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	fprintf(stdout, "I'm number %d of %d MPI processes\n", rank, size);
	fflush(stdout);

	int buf[10];
	memset(buf, 0, sizeof(int) * 10);
	int i;
	MPI_Status status;
	// send and recv
	if (rank % 2 == 0) { // sender
		for (i=0; i<10; i++)
			buf[i] = i*i;
		MPI_Ssend(buf, 10, MPI_INT, ((rank + 1) % size), 999, MPI_COMM_WORLD);
		MPI_Recv(buf, 10, MPI_INT, ((rank + 1) % size), 999, MPI_COMM_WORLD, &status);
		fprintf(stdout, "Proc %d got ", rank);
		for (i=0; i<10; i++)
			fprintf(stdout, "%d ", buf[i]);
		fprintf(stdout, "\n");
		fflush(stdout);
		
	} else { // receiver
		for (i=0; i<10; i++)
			buf[i] = i*i*i;
		MPI_Ssend(buf, 10, MPI_INT, ((rank + size - 1) % size), 999, MPI_COMM_WORLD);
		MPI_Recv(buf, 10, MPI_INT, ((rank + size - 1) % size), 999, MPI_COMM_WORLD, &status);
		fprintf(stdout, "Proc %d got ", rank);
		for (i=0; i<10; i++)
			fprintf(stdout, "%d ", buf[i]);
		fprintf(stdout, "\n");
		fflush(stdout);
	}

	MPI_Finalize();
}
