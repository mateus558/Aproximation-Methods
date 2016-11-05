#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

int main(){
	int N, i, j, e;
	char fname[50], *pos;
	double** matrix = NULL;

	fgets(fname, sizeof(fname), stdin);
	if((pos=strchr(fname, '\n')) != NULL)
    	*pos = '\0';

	if(read_matrix(fname, &matrix, &N) == -1){
		printf("The matrix couldn't be read!\n");
		exit(1);
	}
	
	if(!matrix){
		printf("The matrix couldn't be read!\n");
		exit(1);
	}
	
	for(i = 0; i < N; ++i){
		for(j = 0; j < N+1; ++j){
			printf("%lf ", matrix[i][j]);
		}	
		printf("\n");
	}
	
	return 0;
}

