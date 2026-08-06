#ifndef _UTILS_H_
#define _UTILS_H_
   
typedef double*  VECTOR;
typedef double** MATRIX;

VECTOR Vector( long n );

MATRIX Matrix( long n, long m );

void free_matrix( MATRIX hMatrix );
void free_vector( VECTOR hVector );

#endif /* _UTILS_H_ */
