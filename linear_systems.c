#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EPS 1E-5
#define MAXIT 1E2
#define INF 2000000000
#define NEG_INF -2000000000

int read_matrix(char *fname, double ***matrix, int *N);
double LU_decomposition(double **A, double ***L, double ***U, int N);
int is_singular(double **A, int N);
int relative_error(double* x, double* xk, int N);
double* jacobi_method(double **A, int N, int (*stop_criterion)(double*, double*, int));
double* seidel_method(double **A, int N, int (*stop_criterion)(double*, double*, int));
double* gauss_elimination(double **A, int N);
double* LU_solve(double **A, int N);

int main(){
	int N, M, i, j;
	char fname[50], *pos;
	double **A = NULL;
	double *x = NULL;

	fgets(fname, sizeof(fname), stdin);
	if((pos=strchr(fname, '\n')) != NULL)
    	*pos = '\0';
	
	M = read_matrix(fname, &A, &N);
	
	if(M == -1 || !A){
		printf("The matrix couldn't be read!\n");
		exit(1);
	}
	
	for(i = 0; i < N; ++i){
		for(j = 0; j <= N; ++j){
			printf("%lf ", A[i][j]);
		}	
		printf("\n");
	}
	
	printf("\n");
	//if(is_singular(A, N) == -1) printf("Singular matrix\n");
	x = jacobi_method(A,N, relative_error);
	
	for(i = 0; i < N; ++i)
		printf("%lf ", x[i]);
	printf("\n");

	return 0;
}

int read_matrix(char *fname, double*** matrix, int *N){
	int i, j, lines, cols, size;
	double val;
	FILE *file = NULL;
	
	file = fopen(fname, "r");
	if(!file){
		printf("The file could not be opened!\n");
		return -1;
	}

	fscanf(file, "%d", N);
	lines = (*N);
	cols = lines+1;
	(*matrix) = (double **)malloc(lines * sizeof(double));	
	for(i = 0; i < lines; ++i){ (*matrix)[i] = (double *)malloc(cols * sizeof(double)); }
	
	for(i = 0; i < lines; ++i){
		for(j = 0; j < cols; ++j){
			fscanf(file, "%lf", &(*matrix)[i][j]);
		}	
	}
	
	fclose(file);	
	return cols;
}

double LU_decomposition(double** A, double ***L, double ***U, int N){
	int i, j, k;
	double **Lk = NULL, **Uk = NULL;
	double sumu = 0.0, suml = 0.0, det = 1.0;
	
	Lk = (double **)malloc(N * sizeof(double));	
	for(i = 0; i < N; ++i){ Lk[i] = (double *)malloc(N * sizeof(double)); Lk[i][i] = 1;}
	Uk = (double **)malloc(N * sizeof(double));	
	for(i = 0; i < N; ++i){ Uk[i] = (double *)malloc(N * sizeof(double)); }
	
	for(k = 0; k < N; ++k){
		Uk[k][k] = A[k][k];
		
		for(i = k+1; i < N; ++i){
			Lk[i][k] = A[i][k]/A[k][k];
			Uk[k][i] = A[k][i];
		}
		
		for(i = k+1; i < N; ++i){
			for(j = k+1; j < N; ++j){
				A[i][j] = A[i][j] - Lk[i][k]*Uk[k][j];
			}
		}
	}
	
	*U = Uk;
	*L = Lk;
	
	for(i = 0; i < N; ++i){
		det *= Uk[i][i];
	}

	return det;
}

int is_singular(double **A, int N){
	int i;
	double **L = NULL, **U = NULL;
	
	if(A[0][0] == 0) return -1;
	
	for(i = 1; i <= N; ++i){
		if(LU_decomposition(A, &L, &U, i) == 0)
			return -1;
	}
	
	return 1;
}

int relative_error(double* x, double* xk, int N){
	double diff = 0.0, g_diff = NEG_INF, g = NEG_INF;
	int i;
	
	for(i = 0; i < N; ++i){
		diff = abs(x[i] - xk[i]);
		if(g_diff < diff) g_diff = diff;
		if(g < abs(x[i])) g = abs(x[i]);
	}
	
	return ((g_diff/g) < EPS)?1:0;
}

double* jacobi_method(double **A, int N, int (*stop_criterion)(double*, double*, int)){
	int i, j, k;
	double *x = NULL, *xk = NULL, *b = NULL, *D = NULL, **B = NULL;

	x = malloc(N * sizeof(double));
	for(i = 0; i < N; ++i) x[i] = 0.0;
		
	xk = malloc(N * sizeof(double));
	for(i = 0; i < N; ++i) xk[i] = 0.0;
	
	D = malloc(N * sizeof(double));
	for(i = 0; i < N; ++i) D[i] = A[i][i];
	
	b = malloc(N * sizeof(double));
	for(i = 0; i < N; ++i) b[i] = A[i][N]/D[i];
	
	//Allocate iteration matrix
	B = (double **)malloc(N * sizeof(double));
	//Pre compute iteration matrix
	for(i = 0; i < N; ++i){ 
		B[i] = (double *)malloc(N * sizeof(double)); 
		
		for(j = 0; j < N; ++j){
			if(i != j){
				B[i][j] = -A[i][j]/D[i];
			}else{
				B[i][j] = 0.0;
			}			
		}
	}
	
	k = 0;
	do{
		for(i = 0; i < N; ++i){
			xk[i] = x[i];
			for(j = 0; j < N; ++j){
				x[i] += B[i][j]*x[j];
			}
			x[i] += b[i];
		}
		k++;
	}while(!stop_criterion(x, xk, N));
	
	return x;
}

double* seidel_method(double **A, int N, int (*stop_criterion)(double*, double*, int)){
	int i, j, k;
	double *x = NULL, *xk = NULL, *b = NULL, **L = NULL, **R = NULL;

	x = malloc(N * sizeof(double));
	for(i = 0; i < N; ++i) x[i] = 0.0;
		
	xk = malloc(N * sizeof(double));
	for(i = 0; i < N; ++i) xk[i] = 0.0;
	
	L = malloc(N * sizeof(double));
	R = malloc(N * sizeof(double));
	for(i = 0; i < N; ++i){
		L[i] = (double *)malloc(N * sizeof(double));
	 	R[i] = (double *)malloc(N * sizeof(double));
		
		for(j = 0; j < N; ++j){
			if(i > j){
				L[i][j] = A[i][j];
				R[i][j] = 0.0;
			}else if(i < j){
				R[i][j] = A[i][j];
				L[i][j] = 0.0;
			}
		}
	}
	
	b = malloc(N * sizeof(double));
	for(i = 0; i < N; ++i) b[i] = A[i][N];
	
	k = 0;
	do{
		
		k++;
	}while(!stop_criterion(x, xk, N));	
	
	return x;
}

double* gauss_elimination(double **A, int N){
	int i, j, k;
	double sum = 0.0, *x = NULL, *b = NULL;
	
	x = malloc(N * sizeof(double));
	b = malloc(N * sizeof(double));
	
	for(i = 0; i < N; ++i) b[i] = A[i][N];

	for(i = 0; i < N-1; i++){ 
		for(j = i; j < N; j++){ 
			if(i > j){
			   double ratio = A[j][i]/A[i][i]; 
			   for(k = i; k < N; k++){ 
			        A[j][k] -= (ratio*A[i][k]); 
			        b[j] -= (ratio*b[i]); 
			   } 
			}
       } 
  	}
	
	/*for(i = 0; i < N-1; ++i){
		for(j = i; j < N; ++j){
			double ratio = A[j][i] / A[i][i];
			for(k = i; k < N; ++k){	
				A[j][k] -= A[i][k] * ratio;
				A[i][N] -= ratio;
			}
		}
	}*/
	
	/*for(k = 0; k < N-1; ++k){
		for(i = k + 1; i < N; ++i){
			A[i][N] = A[i][N] - A[k][N] * (A[i][k] / A[k][k]);
			for(j = k; j < N; ++j){
				A[i][j] = A[i][j] - A[k][j] * (A[i][k] / A[k][k]);
			}
			printf("%lf %d %d \n", A[i][k], k, i);
		}
	}*/
	
	for(k = 0; k < N; ++k){
		for(i = 0; i <= N; ++i){
			printf("%lf ", A[k][i]);	
		}
		printf("\n");
	}
	
	for(i = N-1, sum = 0.0; i >= 0; --i){
		x[i] = b[i];
		for(j = i+1; j < N; ++j){
			sum += A[i][j]*x[j];
		}
		x[i] = (b[i] - sum)/A[i][i];
	}
	
	return x;
}

double* LU_solve(double **A, int N){
	int i, j;
	double **L = NULL, **U = NULL, *y = NULL, *x = NULL;
	double det = 1.0, sum;
	
	y = malloc(N * sizeof(double));
	x = malloc(N * sizeof(double));
	
	det = LU_decomposition(A, &L, &U, N);
	
	if(det == 0.0){
		return NULL;
	}
	
	for(i = 0; i < N; ++i){
		sum = 0.0;
		for(j = 0; j < i; ++j){
			sum += L[i][j]*y[j];
		}
		y[i] = A[i][N] - sum;
	}
	
	for(i = N-1; i >= 0; --i){
		sum = 0.0;
		for(j = N-1; j > i; --j){
			sum += U[i][j]*x[j];
		}
		x[i] = (y[i] - sum)/U[i][j];
	}
	
	return x;
}


