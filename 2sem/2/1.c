//calculating sum of 1/i for i = 1 to N
//where N is an argument of prog

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <errno.h>


//#define DEBUG
#ifdef DEBUG
#define DBG
#else
#define DBG if(0)
#endif

const double MAGIC    = 1000;
const int    MAX_PROC = 28;

long   get_N(int, char*[]);
void   close_all_connections(int);
double worker(int, int, int);

int main(int argc, char *argv[]) 
{
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	long N = 0;
	double starttime = 0.0;
	if (rank == 0) 
	{
		N = get_N(argc, argv);
		if (N < 1)
		{
			close_all_connections(size);
			printf("WRONG ARG ABORT");
			return 0;
		}
		starttime = MPI_Wtime();

		//give tasks
		MPI_Request requests     [MAX_PROC];
		MPI_Request requests_recv[MAX_PROC];
		int    flags     [MAX_PROC];
		int    flags_recv[MAX_PROC];
		double res       [MAX_PROC];

		for (int i = 0; i < size - 1; ++i)
		{
			flags        [i] = 0;
			flags_recv   [i] = 0;
			res          [i] = 0;
			requests     [i] = MPI_REQUEST_NULL;
			requests_recv[i] = MPI_REQUEST_NULL;
		}

		for (int i = 1; i < size; ++i)
		{
			MPI_Isend(&N, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &requests[i-1]);
		}
		DBG printf("DBG: all isend\n");

		//collect status and wait for result
		MPI_Status status;	
		int sent   = 0;
		int recved = 0;

		for (int i = 0; i < size - 1; ++i)
		{
			MPI_Irecv(&res[i], 1, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD, &requests_recv[i]);
		}
		DBG printf("DBG: all irecv\n");

		while (sent != size - 1 || recved != size - 1) 
		{
			for (int i = 0; i < size - 1; ++i)
			{
				if (flags[i] == 0)
				{
					MPI_Test(&requests[i], &flags[i], &status);
					if (flags[i] == 1) sent ++;
				}		
				if (flags_recv[i] == 0)
				{
					MPI_Test(&requests_recv[i], &flags_recv[i], &status);
					if (flags_recv[i] == 1) recved ++;
				}
			}
		}

		//calculate result
		double ans = 0;
		for (int i = 0; i < size - 1; ++i) 
		{
			ans += res[i];
			DBG printf("worker #%d: %f\n", i+1, res[i]);
		}
		printf("RESULT: %f\n\n", ans);
	}
	else
	{
		MPI_Status status;
		MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		if (N < 0) {
			MPI_Finalize();
			return -1;
		}	

		double res = worker(N, size, rank);

		MPI_Send(&res, 1, MPI_DOUBLE, 0 , 0, MPI_COMM_WORLD);
	}
	

	if (rank == 0)
	{
		double endtime = MPI_Wtime();
		printf("elapsed time %f\n", (endtime-starttime) * MAGIC);
	}

	MPI_Finalize();
	return 0;
}

double worker(int N, int size, int rank)
{
	DBG printf("DBG: worker %d start with N = %d, size = %d\n", rank, N, size);
	double res = 0;
	int i = rank;

	while (i < N)
	{
		double x = 1.0 / i;
		res += x;
		i += size;
	}
	
	DBG printf("DBG: worker %d end with %f\n", rank, res);
	return res;
} 

long get_N(int argc, char* argv[])
{
        if (argc != 2)
                return -1;

        char* endptr = NULL;
        long ret = strtol(argv[1], &endptr, 10);

        if ((errno == ERANGE && errno == EINVAL) || *endptr != '\0')
        {
                return -1;
        }

        return ret;
}

void close_all_connections(int size)
{
	if (size < 1)
		return;
	
	int errmes = -1;
	for (int i = 1; i < size; ++i) 
	{
		MPI_Send(&errmes, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();
}
