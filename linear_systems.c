#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EPS 1E-5
#define MAXIT 1E2
#define INF 2000000000
#define NEG_INF -2000000000

int read_matrix(char *fname, double ***matrix, int *N);
void print_matrix(double **A, int N);
int is_singular(double **A, int N);
int relative_error(double* x, double* xk, int N);
double relative_error_value(double* x, double* xk, int N);
int line_criterion(double **A, int N);
double** swap_lines(double **A, int N, int L1, int L2);
double** swap_columns(double **A, int N, int L1, int L2);
double LU_decomposition(double **A, double ***L, double ***U, int N);
double* jacobi_method(double **A, int N, int (*stop_criterion)(double*, double*, int));
double* seidel_method(double **A, int N, int (*stop_criterion)(double*, double*, int));
double* gauss_elimination(double **A, int N);
double* back_substitution(double** A, double* b, int N);
double* LU_solve(double **A, int N);

int main(){
	int N, M, i, j;
	char fname[3] = {'i', 'n', '\0'}, *pos;
	double **A = NULL;
	double *x = NULL;

	/*fgets(fname, sizeof(fname), stdin);
	if((pos=strchr(fname, '\n')) != NULL)
    	*pos = '\0';
	*/
	
	M = read_matrix(fname, &A, &N);
	
	if(M == -1 || !A){
		printf("The matrix couldn't be read!\n");
		exit(1);
	}
	
	printf("Original matrix: \n\n");
	
	print_matrix(A, N);
	printf("\n");
	
	x = gauss_elimination(A, N);
	
	if(x){
		printf("\nSolution:\n");
		for(i = 0; i < N; ++i)
			printf("x[%d]: %lf, ", i, x[i]);
		printf("\n");
	}else{
		printf("\nCouldn't find the solution!\n");
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

void print_matrix(double **A, int N){
	int i, j;
	
	for(i = 0; i < N; i++){
    	for(j = 0; j <= N; j++){
    		if(j == N) printf(" | ");
    		printf("%lf ",A[i][j]);
    	}
    	printf("\n");	
   }
}

int is_singular(double **A, int N){
	int i, ret = 1;
	double **L = NULL, **U = NULL;
	
	if(A[0][0] == 0) return -1;
	
	for(i = 1; i <= N; ++i){
		if(LU_decomposition(A, &L, &U, i) == 0){
			ret = 0;
			break;	
		}
	}
	
	for(i = 0; i < N; ++i){
		free(L[i]);
		free(U[i]);
	}
	free(L);
	free(U);
	
	return ret;
}

int line_criterion(double **A, int N){
	int i, j, sign = 1;
	double alpha = 0.0, sum = 0.0, max_alpha = NEG_INF, ratio;
	
	for(i = 0; i < N; ++i){
		for(j = 0; j < N; ++j){
			if(i != j){
				ratio = A[i][j] / A[i][i];
				sign = (ratio < 0)?-1:1;
				alpha += ratio * sign;
			}
		}
		if(alpha > max_alpha) max_alpha = alpha;
	}
	
	return (max_alpha < 1)?1:-1;
}

double relative_error_value(double* x, double* xk, int N){
	double g_diff = NEG_INF, g = NEG_INF;
	double *diff = NULL;
	int i = 0;
	
	diff = malloc(N * sizeof(double));
	
	for(i = 0; i < N; ++i){
		diff[i] = fabs(x[i] - xk[i]);
	}
	for(i = 0; i < N; ++i){
		if(diff[i] > g_diff) g_diff = diff[i];
		if(x[i] > g) g = x[i];
	}

	free(diff);

	return fabs(g_diff/g);
}

int relative_error(double* x, double* xk, int N){
	double g_diff = NEG_INF, g = NEG_INF;
	double *diff = NULL;
	int i = 0;
	
	diff = malloc(N * sizeof(double));
	
	for(i = 0; i < N; ++i){
		diff[i] = fabs(x[i] - xk[i]);
	}
	for(i = 0; i < N; ++i){
		if(diff[i] > g_diff) g_diff = diff[i];
		if(x[i] > g) g = x[i];
	}

	free(diff);

	return ((fabs(g_diff/g)) < EPS)?1:0;
}

double** swap_lines(double **A, int N, int L1, int L2){
	int i, j, k;
	int **P = NULL;
	double sum = 0.0, temp, **res = NULL;
	
	res = (double **)malloc(N * sizeof(double));
	for(i = 0; i < N; ++i){ res[i] = (double *)malloc(N+1 * sizeof(double)); }
	P = (int **)malloc(N * sizeof(int));
	for(i = 0; i < N; ++i){ P[i] = (int *)malloc(N * sizeof(int)); }
	for(i = 0; i < N; ++i){
		for(j = 0; j < i; ++j){
			P[i][j] = 0;
		}
		res[i][N] = A[i][N];
		P[i][i] = 1;
		for(j = i + 1; j < N; ++j){
			P[i][j] = 0;
		}
	}
	
	P[L1][L1] = 0;
	P[L1][L2] = 1;
	P[L2][L2] = 0;
	P[L2][L1] = 1;
	
	for(i = 0; i < N; i++){
    	for(j = 0; j < N; j++){
    		for(k = 0; k < N; k++){
    			sum += P[i][k] * A[k][j];
        	}
 
        	res[i][j] = sum;
        	sum = 0.0;
    	}
	}
	
	res[L1][N] = A[L2][N];
	res[L2][N] = A[L1][N];
	
	return res;
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

double* jacobi_method(double **A, int N, int (*stop_criterion)(double*, double*, int)){
	int i, j, k;
	double *x = NULL, *xk = NULL, *b = NULL;
	double s;
	
	x = malloc(N * sizeof(double));
	for(i = 0; i < N; ++i) x[i] = 0.0;
	
	if(line_criterion(A, N) == -1) return NULL;
		
	xk = malloc(N * sizeof(double));
	for(i = 0; i < N; ++i) xk[i] = 0.0;
		
	b = malloc(N * sizeof(double));
	for(i = 0; i < N; ++i) b[i] = A[i][N];
	
	printf("\n");
	for(i = 0; i <= N; ++i)	printf("-----------------");
	printf("\n");
	printf(" itr	");
	for(i = 0; i < N; ++i){
		printf("x[%d]		", i);
	}
	printf("  error");
	printf("\n");
	for(i = 0; i <= N; ++i)	printf("-----------------");
	printf("\n");	
	
	k = 0;
	
	printf(" %d	", k);
	for(i = 0; i < N; ++i){
		printf("%lf	", x[i]);
	}
	printf("    -");
	printf("\n");
	
	do{
		printf(" %d	", k+1);
		s = 0.0;
		for(i = 0; i < N; ++i)
			xk[i] = x[i];
	
		for(i = 0; i < N; ++i){
			for(j = 0; j <= i - 1; ++j){
				s += A[i][j] * xk[j];
			}
			for(j = i + 1; j < N; ++j){
				s += A[i][j] * xk[j];
			}
			
			x[i] = (b[i] - s)/A[i][i];
			printf("%lf	", x[i]);
		}
		k++;
		printf("%lf", relative_error_value(x, xk, N));
		printf("\n");
	}while(!stop_criterion(x, xk, N) && k < MAXIT);

	free(xk);
	free(b);

	return x;
}

double* seidel_method(double **A, int N, int (*stop_criterion)(double*, double*, int)){
	int i, j, k;
	double *x = NULL, *xk = NULL, *b = NULL;
	double s;
	
	x = malloc(N * sizeof(double));
	for(i = 0; i < N; ++i) x[i] = 0.0;
	
	if(line_criterion(A, N) == -1) return NULL;
	
	xk = malloc(N * sizeof(double));
	for(i = 0; i < N; ++i) xk[i] = 0.0;
		
	b = malloc(N * sizeof(double));
	for(i = 0; i < N; ++i) b[i] = A[i][N];
	
	printf("\n");
	for(i = 0; i <= N; ++i)	printf("-----------------");
	printf("\n");
	printf(" itr	");
	for(i = 0; i < N; ++i){
		printf("x[%d]		", i);
	}
	printf("  error");
	printf("\n");
	for(i = 0; i <= N; ++i)	printf("-----------------");
	printf("\n");	

	k = 0;
	
	printf(" %d	", k);
	for(i = 0; i < N; ++i){
		printf("%lf	", x[i]);
	}
	printf("    -");
	printf("\n");	
	
	do{
		s = 0.0;
		printf(" %d	", k+1);
		for(i = 0; i < N; ++i)
			xk[i] = x[i];
			
		for(i = 0; i < N; ++i){
			
			for(j = 0; j <= i - 1; ++j){
				s += A[i][j] * x[j];
			}
			for(j = i + 1; j < N; ++j){
				s += A[i][j] * xk[j];
			}
			
			x[i] = (b[i] - s)/A[i][i];
			printf("%lf	", x[i]);
		}
		printf("%lf", relative_error_value(x, xk, N));
		printf("\n");
		k++;
	}while(!stop_criterion(x, xk, N) && k < MAXIT);	
	
	free(xk);
	free(b);
	
	return x;
}

double* back_substitution(double** A, double* b, int N){
	int i, j;
	double sum = 0.0, *x = NULL;
	
	x = malloc(N * sizeof(double));
	
	x[N-1] = b[N-1]/A[N-1][N-1];
	
	for(i = N-1, sum = 0.0; i >= 0; --i){
		sum = b[i];
		for(j = i+1; j < N; ++j){
			sum -= A[i][j]*x[j];
		}
		x[i] = sum/A[i][i];
	}
	
	return x;
}

double* gauss_elimination(double **A, int N){
	int i, j, k, l = -1, flag = 1;
	double ratio, sum = 0.0, *x = NULL, *b = NULL;
	
	x = malloc(N * sizeof(double));
	b = malloc(N * sizeof(double));
	
	for(i = 0; i < N; ++i){
		for(j = 0; j < i; ++j){
			if(fabs(A[j][i]) > fabs(A[i][i]) && A[j][i] != 0) l = j;
		}
		for(j = i+1; j < N; ++j){
			if(fabs(A[j][i]) > fabs(A[i][i]) && A[j][i] != 0) l = j;
		}

		if(l != -1)
			A = swap_lines(A, N, i, l);
		l = -1;
	}
	
	for(i = 0; i < N; ++i) b[i] = A[i][N];
	
	printf("\nPivoted Matrix: \n");
	print_matrix(A, N);
	
	for(k = 0; k < N - 1; ++k){
		for(i = k + 1; i < N; ++i){
			if(A[k][k] == 0){
				for(j = 0; j < k; ++j){
					if(A[j][k] != 0 && fabs(A[j][k]) > fabs(A[l][k])) l = j;
				}
				
				for(j = k+1; j < N; ++j){
					if(A[j][k] != 0 && fabs(A[j][k]) > fabs(A[l][k])) l = j;
				}
				
				printf("%d %d\n", k, l);
				if(l != -1)
					A = swap_lines(A, N, k, l);
				l = -1;
			}
			ratio = A[i][k]/A[k][k];
			A[i][k] -= A[i][k];
			for(j = k + 1; j < N; ++j){
				A[i][j] -= ratio * A[k][j];
			}
			b[i] -= ratio * b[k];
		}
	}
	
	printf("\nScalonated matrix:\n");
	print_matrix(A, N);
	
	for(i = 0; i < N; ++i){
		if(A[N-1][i] != 0){ flag = 0; break; }
	}	
	
	if(flag){
		printf("\nSingular matrix.\n");
		return NULL;
	}
	
	x = back_substitution(A, b, N);
	
	free(b);
	
	return x;
}

double* LU_solve(double **A, int N){
	int i, j;
	double **L = NULL, **U = NULL, *y = NULL, *x = NULL;
	double det = 1.0, sum = 0.0;
	
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
	
	for(i = 0; i < N; ++i){
		free(U[i]);
		free(L[i]);
	}
	free(U);
	free(L);
	free(y);
	
	return x;
}


