1.o: 1.c
	mpicc -Wall -Wextra -std=c99 -lm 1.c -o 1.o
	qsub job1
test:
	rm 81410.e*
	qsub job1
clean:
	# rm *.o
	rm 81410.*
	rm core.*
	rm *.swp
