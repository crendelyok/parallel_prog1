//cirlce transmition
#include <stdio.h>
#include <mpi.h>

//const double MAGIC = 0;

int main(int argc, char *argv[]) 
{
	int size, rank;
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int a = 0;
	//double starttime = 0.0;
	if (rank == 0) 
	{
		//double starttime = MPI_Wtime();
		//printf("#%d: send int = %d\n\n", rank, a);
		//fflush(stdout);

		MPI_Send(&a, 1, MPI_INT, 1 % size, 0, MPI_COMM_WORLD);

		MPI_Status status;
		MPI_Recv(&a, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, &status);
		//printf("#%d: recv int = %d\n", rank, a);
	}
	else
	{
		MPI_Status status;
		MPI_Recv(&a, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);
		//printf("#%d: recv int = %d\n", rank, a);

		a++;		
		//printf("#%d: int++\n", rank);

		//printf("#%d: sent int = %d\n\n", rank, a);
		//fflush(stdout);
		MPI_Send(&a, 1, MPI_INT, (rank + 1) % size , 0, MPI_COMM_WORLD);
	}
	

	if (rank == 0)
	{
		//double endtime = MPI_Wtime();
		printf("END OF CIRCLE, int = %d\n", a);
		//printf("elapsed time %f\n", (endtime-starttime));
	}

	MPI_Finalize();
	return 0;
}
