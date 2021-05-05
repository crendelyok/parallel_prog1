//lab1

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define DEBUG
#ifdef DEBUG
#define DBG
#else
#define DBG if(0)
#endif


const double PI = 3.141592653589;
const double MAGIC    = 1000;
const int    MAX_PROC = 28;
const int    T_STEPS  = 2000;
const int    X_STEPS  = 2000;
const double MAX_TIME = 1.0; // t is from 0 to MAX_TIME
const double MAX_X    = 1.0; // x is from 0 to MAX_X

//long   get_N(int, char*[]);
void   close_all_connections(int);
double worker(int, int, int);
double FuncU_T0(double);
double FuncU_X0(double);
double FuncF(double, double);
double FuncA(); //let make it = const

int main(int argc, char *argv[]) 
{
	
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	//assigning src and dest ranks
	int dest = rank + 1;
	if (rank == size-1)
		dest = MPI_PROC_NULL;
	int src = rank - 1;
	if (src == -1)
		src = MPI_PROC_NULL;

	//DBG printf("size %d, rank %d, dest %d, src %d\n\n", size, rank, dest, src);

	double t_step = MAX_TIME / T_STEPS;
	double x_step = MAX_X / X_STEPS;

	if (rank == 0) 
	{
		double starttime = MPI_Wtime();

		//giving tasks
		MPI_Request requests     [MAX_PROC];
		//MPI_Request requests_recv[MAX_PROC];
		//int flags     [MAX_PROC];
		//int flags_recv[MAX_PROC];

		for (int i = 0; i < size - 1; ++i)
		{
			//flags        [i] = 0;
			//flags_recv   [i] = 0;
			//res          [i] = 0;
			requests     [i] = MPI_REQUEST_NULL;
			//requests_recv[i] = MPI_REQUEST_NULL;
		}

		int worker_starts[MAX_PROC]; // left step-boundary x-axis 
		worker_starts[0] = 0;
		/* right - left step-boundaries x-axis aka num of operations per worker */
		int worker_sizes [MAX_PROC];
		int iteration = X_STEPS;
		for (int i = 0; i < size; ++i)
		{
			/*
			worker_starts[i] = X_STEPS / size * i;
			worker_sizes [i] = X_STEPS / size;
			
			if (i == size - 1)
				worker_sizes[i] = X_STEPS - worker_starts[i];

			*/

			if (iteration % (size - i) != 0)
				worker_sizes[i] = iteration / size + 2;
			else 
				worker_sizes[i] = iteration / (size - i) + 1;
			iteration -= worker_sizes[i] - 1;
			DBG printf("w%d: start=%d, size=%d\n", i, worker_starts[i], worker_sizes[i]);
			if (i != 0) {
				MPI_Isend(&worker_starts[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &requests[i-1]);
				MPI_Isend(&worker_sizes [i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &requests[i-1]);
			}

			if (i != size - 1)
				worker_starts[i+1] = worker_starts[i] + worker_sizes[i] - 1;
		}
		DBG printf("DBG: master gave tasks\n");


		// master act like a worker
		

		MPI_Request request; // empty request
		MPI_Status status;   // empty status
		double** result_table = (double**) calloc((X_STEPS + 1), sizeof(double*));
		for (int i = 0; i <= X_STEPS; ++i) 
			result_table[i] = (double*) calloc((T_STEPS + 1), sizeof(double));

		for (int i = 0; i<= T_STEPS; ++i) 
			result_table[0][i] = FuncU_T0(i * t_step);	
		for (int i = 0; i <= X_STEPS; ++i)
			result_table[i][0] = FuncU_X0(i * x_step);

		for (int t = 0; t < T_STEPS; ++t)
		{
			for (int x = 1; x < worker_sizes[0]; ++x)
				result_table[x][t+1] = ((2 * x_step * t_step) / (x_step + t_step * FuncA())) *
						       (FuncF((x + 0.5) * x_step, (t + 0.5) * t_step) +
						       FuncA() * (result_table[x-1][t+1] - result_table[x][t] +
					   	       result_table[x-1][t]) / (2 * x_step) + 
						       (result_table[x-1][t] - result_table[x-1][t+1] + result_table[x][t]) /
						       (2 * t_step));
			MPI_Isend(&result_table[worker_sizes[0]-1][t+1], 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &request);
			//DBG printf("master isent to %d: %.3f\n", dest, result_table[worker_sizes[0]-1][t+1]);
		}
		DBG printf("DBG: master did his job B-)\n\n");

		/*
		if (size == 1) {
			for (int x = 0; x <= X_STEPS; ++x) 
			{
				for (int t = 0; t <= T_STEPS; ++t) 
					DBG printf("%.3f ", result_table[x][t]);
				DBG printf("\n");
			}
		}		
		*/

		// collect status and wait for result

		for (int i = 0; i < size - 1; ++i)
			for (int j = worker_starts[i+1] + 1; j < worker_starts[i+1] + worker_sizes[i+1]; ++j) {
				//DBG printf("master final waits from %d\n", i);
				//printf("pointer %p\n", result_table[j]);
				MPI_Recv(result_table[j], T_STEPS+1, MPI_DOUBLE, i + 1, 3, MPI_COMM_WORLD, &status);
				//DBG printf("master final recved from %d\n", i+1);
			}
		double endtime = MPI_Wtime();
		printf("elapsed time: %f\n\n", endtime-starttime);

		// free allocated memory 
		for (int i = 0; i <= X_STEPS; ++i) 
			free(result_table[i]);
		free(result_table);
	}
	else
	{
		// worker
		MPI_Status status;
		MPI_Request request;
		int this_start, this_size;
		MPI_Recv(&this_start, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&this_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		
		//DBG printf("DBG: worker%d recv: start%d, size%d\n", rank, this_start, this_size);
		//DBG printf("DBG: w#%d src: %d, dest %d\n", rank, src, dest);

		double** table = (double**) calloc(this_size, sizeof(double*));
		for (int i = 0; i < this_size; i++)
			table[i] = (double*) calloc(T_STEPS + 1, sizeof(double));

		for (int i = 0; i < this_size; ++i)
			table[i][0] = FuncU_X0(this_start * x_step);

		for (int t = 0; t < T_STEPS; ++t)
		{
			//DBG printf("w#%d waits from %d\n", rank, src);
			MPI_Recv(&table[0][t+1], 1, MPI_DOUBLE, src, 0, MPI_COMM_WORLD, &status);
			//DBG printf("w#%d recved from %d\n", rank, src);
			for (int x = 1; x < this_size; ++x) 
			{
				int abs_x = x + this_start; // absolute x (as it will be in result_table)
				table[x][t+1] = ((2 * x_step * t_step) / (x_step * t_step * FuncA())) *
  					        (FuncF((abs_x + 0.5) * x_step, (t + 0.5) * t_step) +
						FuncA() * (table[x-1][t+1] - table[x][t] + table[x-1][t]) / 
						(2 * x_step) + (table[x-1][t] - table[x-1][t+1] + table[x][t]) /
						(2 * x_step));
			}

			MPI_Isend(&table[this_size-1][t+1], 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &request);
			//DBG printf("w#%d sent to %d\n", rank, dest);
		}
			
		for (int x = 1; x < this_size; ++x) {
			MPI_Send(table[x], T_STEPS+1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
			//DBG printf("w#%d final sent to master\n", rank);
		}

		for (int i = 0; i < this_size; ++i)
			free(table[i]);
		free(table);

		DBG printf("DBG: worker%d finished\n", rank);
	}
	
	MPI_Finalize();
	return 0;
}

double FuncU_T0(double t)
{
	return exp(-t);
}

double FuncU_X0(double x)
{
	return cos(PI * x);
}

double FuncF(double x, double t)
{
	return x + t;
}

double FuncA()
{
	return 2;
}

/*
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
*/
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
