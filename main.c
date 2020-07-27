#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "matrix.h"

const size_t side = 1000;

typedef void (*mm_mult_t)(matrix_t result,
			  const matrix_t left_source,
			  const matrix_t right_source);
void time_mm_mult(mm_mult_t multiplication,
		  matrix_t result,
		  matrix_t left_source,
		  matrix_t right_source,
		  size_t num_tests);

int main()
{
	srand48(time(0));
	matrix_t A = new_random_matrix(side);
	matrix_t B = new_random_matrix(side);
	matrix_t C = new_random_matrix(side);
	printf("Naive matrix matrix multiplication\n");
	time_mm_mult(mm_mult_naive,C,A,B,10);
	printf("Blocked matrix matrix multiplication\n");
	time_mm_mult(mm_mult_blocked,C,A,B,10);
	printf("Row major column major matrix multiplication\n");
	time_mm_mult(mm_mult_rm_cm,C,A,B,10);
	printf("Blocked row major column major matrix matrix multiplication\n");
	time_mm_mult(mm_mult_blocked_rm_cm,C,A,B,10);
	printf("Blocked vectorized row major column major matrix matrix multiplication\n");
	time_mm_mult(mm_mult_blocked_rm_cm_vec,C,A,B,10);
	printf("Blocked parallel vectorized row major column major matrix matrix multiplication\n");
	time_mm_mult(mm_mult_blocked_rm_cm_vec_parallel,C,A,B,10);
	free_matrix(A);
	free_matrix(B);
	free_matrix(C);
	return EXIT_SUCCESS;
}

void time_mm_mult(mm_mult_t multiplication,
		  matrix_t result,
		  matrix_t left_source,
		  matrix_t right_source,
		  size_t num_tests)
{
	for (size_t i = 0; i<num_tests; i++)
	{
		printf("Test %lu:",i+1);
		fflush(stdout);
		struct timespec t_start,t_end;
		clock_gettime(CLOCK_REALTIME,&t_start);
		multiplication(result,left_source,right_source);
		clock_gettime(CLOCK_REALTIME,&t_end);
		double elapsed_time = 
			(t_end.tv_sec-t_start.tv_sec)*1e6+
			(t_end.tv_nsec-t_start.tv_nsec)*1e-3;
		printf(" %10.3lf µs\n",elapsed_time);
	}
}
