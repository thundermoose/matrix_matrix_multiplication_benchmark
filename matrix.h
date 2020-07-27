#ifndef __MATRIX__
#define __MATRIX__

#include <stdlib.h>

struct _matrix_;
typedef struct _matrix_ *matrix_t;

matrix_t new_random_matrix(size_t side);

void print_matrix(const matrix_t matrix);

void mm_mult_naive(matrix_t result,
		   const matrix_t left_source,
		   const matrix_t right_source);

void mm_mult_blocked(matrix_t result,
		     const matrix_t left_source,
		     const matrix_t right_source);

void mm_mult_rm_cm(matrix_t result,
		   const matrix_t left_source,
		   const matrix_t right_source);

void mm_mult_blocked_rm_cm(matrix_t result,
			   const matrix_t left_source,
			   const matrix_t right_source);
void mm_mult_blocked_rm_cm_vec(matrix_t result,
			       const matrix_t left_source,
			       const matrix_t right_source);
void mm_mult_blocked_rm_cm_vec_parallel(matrix_t result,
					const matrix_t left_source,
					const matrix_t right_source);

void free_matrix(matrix_t matrix);

#endif
