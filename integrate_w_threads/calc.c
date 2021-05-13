#define _GNU_SOURCE

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/sysinfo.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <sched.h>

#define DEBUG 

#ifdef DEBUG
#define DBG
#else
#define DBG if(0)
#endif

int validatearg(char*);
int get_n_proc();
void* thread_routine(void*);
//void* fiction_routine(void*);
double get_dx(int);
double func_dfdx2(double);
double max_f_dfdx2(double);
double get_extrema(long);
void print_with_precision(double, int);
int count_zeros();


const double PI = 3.1415926535897932384626433832795028841;
const int MAX_PRECISION = 38;
const double MAX_DFDX = 1000000000000000.0; // max df/dx over [0.001, 1] ~ e+15
//const double MAX_DFDX = 100000000000000000000.0; // max df/dx over [0.0001, 1] ~ e+20


const int MAX_PROC = 100;
//double func_dx = 0.0001;
const double func_l = 0.001;
const double func_r = 1;
double segment = 0;
double e = 0;
int zeros; // number of zeros in [left; right]


typedef struct routine_arg {
	double left;
	double right;
	double max_x; //maxima of x over [left, right]
	double ans;
	int tid;
}arg_t;

arg_t* routine_arg_init(int);

//max precision = 38 because of PI
int main(int argc, char* argv[]) {

	if (argc != 3) {
		printf("I need 2 args\n");
		return -1;
	}

	int n_threads = validatearg(argv[1]);
	int precision = validatearg(argv[2]);
	if (precision > MAX_PRECISION) precision = MAX_PRECISION;

	int n_proc = get_nprocs();
	
	e = get_dx(precision);
	e/= n_threads * 10;
	DBG printf("epsilon: %f\n", e);


	DBG printf("Will use %d threads\n", n_threads);

	//int n_fict_thread = (n_threads > n_proc) ? n_threads : n_proc;
	
	pthread_t* tids = calloc(n_threads, sizeof(pthread_t));
	if (tids == NULL) {
		perror("calloc\n");
		exit(1);
	}

	pthread_attr_t attr;
	pthread_attr_init(&attr);

	//segment = (func_r - func_l) / n_threads;	
	arg_t* args = routine_arg_init(n_threads);

	cpu_set_t cpu_set;

	for (int i = 0; i < n_threads; ++i) {
		CPU_ZERO(&cpu_set);
		CPU_SET(i % n_proc, &cpu_set); 

		int ret_code = 0;
		//0 on success
		ret_code = pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t),
				&cpu_set);
		
		if (ret_code != 0) {
			perror("pthread_attr_setaffinity\n");
			exit(1);
		}

		ret_code = pthread_create(&tids[i], &attr, 
				thread_routine, &args[i]);

		DBG printf("created thread#%d\n", i);

		if (ret_code != 0) {
			perror("pthread_create\n");
			exit(1);
		}
	}

	for (int i = 0; i < n_threads; ++i) {
		int ret_code = 0;
		ret_code = pthread_join(tids[i], NULL);
		if (ret_code != 0) {
			perror("pthread_join\n");
			exit(1);
		}
		//DBG printf("thread %d joined\n", i);
	}

	//calc ans
	double ans = 0;
	for (int i = 0; i < n_threads; ++i) {
		ans += args[i].ans;
		}

	print_with_precision(ans, precision);

	pthread_attr_destroy(&attr);

	free(tids);
	free(args);

	return 0;
}

arg_t* routine_arg_init(int n_threads) {
	arg_t* args = (arg_t*) calloc(n_threads, sizeof(arg_t));
	
	DBG printf("number of zeros %d\n", zeros);

	double x = func_l;
	//debug case

	long cur_max_k = 0;
	double cur_max_x = get_extrema(cur_max_k);
	DBG printf("first extrema %.3lf\n", cur_max_x); 
	for (int i = 0; i < n_threads; ++i) {
		args[i].tid = i;
		args[i].left = x;
		if (x * 10 >= func_r) {
			args[i].right = (func_r - args[i].left) / 2;
		}
		else {
			args[i].right = x * 10;
		}
			x = args[i].right;

		if (i == n_threads-1)
			args[i].right = func_r;
		while (get_extrema(cur_max_k) >= args[i].right)
			cur_max_k++;
		while (get_extrema(cur_max_k) >= args[i].left)
			cur_max_k++;
	
		args[i].max_x = get_extrema(cur_max_k - 1);

		DBG printf("%2d: [%.8f, %.8f] with extrema %.5f\n", 
			i, args[i].left, args[i].right, args[i].max_x);

	}

	return args;
}

double func(double x) {
	return sin(1 / x) / x;
}

double func_dfdx2(double x) {
	return abs(((4*x*cos(1/x) + (-1 + 2*x*x)*sin(1/x)) / x*x*x*x*x));
}

double get_extrema(long k) {
	return (2 / (PI * (4 * k + 1)));
}

double max_f_dfdx2(double x) {
	return 5/(x*x*x*x*x);
}

void* thread_routine(void* param) {

	double l = ((arg_t*)param) -> left;
	double r = ((arg_t*)param) -> right;
	//double max_x = ((arg_t*)param) -> max_x;
	double max_dfdx2 = max_f_dfdx2(l);
	double module = r - l;
	long n_dx = (long) round(sqrt(max_dfdx2 * module * module * module / 12 / e));

	double ans = 0.0;
	int tid = ((arg_t*)param) -> tid;
	
	if (n_dx <= 1) {
		DBG printf("WARNING: N_DX <= 1!!\n");
		n_dx = 1;
	}

	double dx = (r-l)/n_dx;

	DBG printf("tid:%d ,n_dx:%ld, max of df/fx^2 is %.1lf\n", tid, n_dx, max_dfdx2);

	double func_left = func(l);
	double func_right = 0.0;
	for (long i = 1; i < n_dx; i++) {
		func_right = func(l + dx * i);
		ans += (func_right + func_left) * dx / 2;
		func_left = func_right;
	}
	// 

	((arg_t*)param) -> ans = ans;

	DBG printf("worker#%d:", tid);
	DBG printf("%f\n", ans);
	//DBG print_with_precision(ans, MAX_PRECISION);
	pthread_exit(NULL);
}

int validatearg(char* argv) {
	char* endptr = NULL;
	int val = strtol(argv, &endptr, 10);

	if (errno == ERANGE || errno == EINVAL || *endptr != '\0') {
		perror("strtol");
		exit(-1);
	}

	if (val < 1) {
		printf("val < 1");
		exit(-1);
	}
	return val;	
}

double get_dx(int precision) {
	double ret = 1.0;
	while (precision > 0) {
		ret /= 10;
		precision --;
	}
	return ret;
}

void print_with_precision(double ans, int precision) {
	printf("%ld.", (long)ans);
	if (ans < 0) ans = -ans;
	ans -= (long)ans;
	for (;precision > 0; precision--) { 
		ans *= 10;
		printf("%ld", ((long)ans) % 10);
		ans -= (long)ans;
	}
	printf("\n");
	return;
}

int count_zeros() {
	int i = 1;
	while(1 / PI / i >= func_l && 1 /PI / i <= func_r) ++i;
	return i-1;
}
