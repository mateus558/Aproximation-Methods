#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EPS 1E-5
#define MAXIT 1E2

int read_matrix(char *fname, double ***matrix, int *N);
int LU_decomposition(double** A, double ***L, double ***U, int N);

int main(){
	int N, M, i, j, e;
	char fname[50], *pos;
	double **matrix = NULL, **L = NULL, **U = NULL;
	double detA = 1.0;

	fgets(fname, sizeof(fname), stdin);
	if((pos=strchr(fname, '\n')) != NULL)
    	*pos = '\0';
	
	M = read_matrix(fname, &matrix, &N);
	if(M == -1){
		printf("The matrix couldn't be read!\n");
		exit(1);
	}
	
	if(!matrix){
		printf("The matrix couldn't be read!\n");
		exit(1);
	}
	
	for(i = 0; i < N; ++i){
		for(j = 0; j < M; ++j){
			printf("%lf ", matrix[i][j]);
		}	
		printf("\n");
	}
	
	printf("\n");
	
	LU_decomposition(matrix, &L, &U, N);
	
	for(i = 0; i < N; ++i){
		for(j = 0; j < N; ++j){
			printf("%lf ", U[i][j]);
		}	
		detA *= U[i][i];
		printf("\n");
	}
	
	if(detA == 0){
		printf("Singular matrix.\n");
		exit(1);
	}
	
	
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

int LU_decomposition(double** A, double ***L, double ***U, int N){
	int i, j, k;
	double **Lk = NULL, **Uk = NULL;
	double sumu = 0.0, suml = 0.0;
	
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
	
	return 1;
}

