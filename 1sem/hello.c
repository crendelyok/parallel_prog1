#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) 
{
	int commsize, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &commsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	printf("task number = %d; rank = %d \n", commsize, rank);
	MPI_Finalize();
	return 0;
}
