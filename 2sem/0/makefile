0.o: 0.c
	mpicc -Wall -Wextra -std=c99 0.c -o 0.o
	qsub job0
test:
	rm 81410.e*
	qsub job1
clean:
	rm *.o
	rm 81410.*
	rm core.*
