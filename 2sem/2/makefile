exp.o: exp.c
	mpicc -Wall -Wextra -std=c99 exp.c -o exp.o
	qsub job_exp
test:
	rm 81410.e*
	qsub job_exp
clean:
	rm *.o
	rm 81410.*
	rm core.*
