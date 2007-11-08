#include <R.h>
#include <Rmath.h>
#include "my_cholesky.h"

/* ********** MODIFIED CHOLESKY FROM THE PACKAGE SURVIVAL ************ */

/* Cholesky factorisation and linear system solver:
	modified versions of functions from the package survival by Terry Therneau
	modification is that I don't use ragged arrays */

/*  SCCS @(#)cholesky2.c	5.2 10/27/98
** subroutine to do Cholesky decompostion on a matrix: C = FDF'
**   where F is lower triangular with 1's on the diagonal, and D is diagonal
**
** arguments are:
**     n         the size of the matrix to be factored
**     **matrix  a ragged array containing an n by n submatrix to be factored
**     toler     the threshold value for detecting "singularity"
**
**  The factorization is returned in the lower triangle, D occupies the
**    diagonal and the upper triangle is left undisturbed.
**    The lower triangle need not be filled in at the start.
**
**  Return value:  the rank of the matrix (non-negative definite), or -rank
**     it not SPD or NND
**
**  If a column is deemed to be redundant, then that diagonal is set to zero.
**
**   Terry Therneau
*/

int my_cholesky2(double *matrix, int *n, double *toler)
{
	double temp;
	int  i,j,k;
	double eps, pivot;
	int rank;
	int nonneg;

	nonneg=1;
	eps =0;
	for (i=0; i<*n; i++) {
		if (matrix[i+i**n] > eps)  eps = matrix[i+i**n];
		for (j=(i+1); j<*n; j++)  matrix[j+i**n] = matrix[i+j**n];
	}
	eps *= *toler;

	rank =0;
	for (i=0; i<*n; i++) {
		pivot = matrix[i+i**n];
		if (pivot < eps) {
			matrix[i+i**n] =0;
			if (pivot < -8*eps) nonneg= -1;
		}
		else  {
			rank++;
			for (j=(i+1); j<*n; j++) {
				temp = matrix[j+i**n]/pivot;
				matrix[j+i**n] = temp;
				matrix[j+j**n] -= temp*temp*pivot;
				for (k=(j+1); k<*n; k++) matrix[k+j**n] -= temp*matrix[k+i**n];
			}
		}
	}
	return(rank * nonneg);
}

/*  SCCS @(#)chsolve2.c	5.2 10/27/98
** Solve the equation Ab = y, where the cholesky decomposition of A and y
**   are the inputs.
**
** Input  **matrix, which contains the chol decomp of an n by n
**   matrix in its lower triangle.
**        y[n] contains the right hand side
**
**  y is overwriten with b
**
**  Terry Therneau
*/

void my_chsolve2(double *matrix, int *n, double *y)
{
	register int i,j;
	register double temp;

     /*
	** solve Fb =y
     */
	for (i=0; i<*n; i++) {
		temp = y[i] ;
		for (j=0; j<i; j++)
			temp -= y[j] * matrix[i+j**n] ;
		y[i] = temp ;
	}
     /*
	** solve DF'z =b
     */
	for (i=(*n-1); i>=0; i--) {
		if (matrix[i+i**n]==0)  y[i] =0;
		else {
			temp = y[i]/matrix[i+i**n];
			for (j= i+1; j<*n; j++)
				temp -= y[j]*matrix[j+i**n];
			y[i] = temp;
		}
	}
}

/*  SCCS @(#)chinv2.c	5.3 07/15/99
** matrix inversion, given the FDF' cholesky decomposition
**
** input  **matrix, which contains the chol decomp of an n by n
**   matrix in its lower triangle.
**
** returned: the upper triangle + diagonal contain (FDF')^{-1}
**            below the diagonal will be F inverse
**
**  Terry Therneau
*/

void my_chinv2(double *matrix, int *n)
{
	register double temp;
	register int i,j,k;

     /*
	** invert the cholesky in the lower triangle
	**   take full advantage of the cholesky's diagonal of 1's
     */
	for (i=0; i<*n; i++){
		if (matrix[i+i**n] >0) {
			matrix[i+i**n] = 1/matrix[i+i**n];   /*this line inverts D */
			for (j= (i+1); j<*n; j++) {
				matrix[j+i**n] = -matrix[j+i**n];
				for (k=0; k<i; k++)     /*sweep operator */
					matrix[j+k**n] += matrix[j+i**n]*matrix[i+k**n];
			}
		}
	}

     /*
	** lower triangle now contains inverse of cholesky
	** calculate F'DF (inverse of cholesky decomp process) to get inverse
	**   of original matrix
     */
	for (i=0; i<*n; i++) {
		if (matrix[i+i**n]==0) {  /* singular row */
			for (j=0; j<i; j++) matrix[j+i**n]=0;
			for (j=i; j<*n; j++) matrix[i+j**n]=0;
		}
		else {
			for (j=(i+1); j<*n; j++) {
				temp = matrix[j+i**n]*matrix[j+j**n];
				if (j!=i) matrix[i+j**n] = temp;
				for (k=i; k<j; k++)
					matrix[i+k**n] += temp*matrix[j+k**n];
			}
		}
	}
}

