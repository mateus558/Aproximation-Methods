#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double** read_matrix(char *fname, double** matrix1, int *N){
	int i, j, e, size;
	double **matrix = NULL, val;
	char *pos;
	FILE *file = NULL;
	
	fgets(fname, sizeof(fname), stdin);
	if((pos=strchr(fname, '\n')) != NULL)
    	*pos = '\0';
	
	file = fopen(fname, "r");
	if(!file){
		printf("The file could not be opened!\n");
		return NULL;
	}

	e = fscanf(file, "%d", N);
	size = (*N);
	matrix = (double **)malloc(size * sizeof(double));	
	for(i = 0; i < size; ++i){ matrix[i] = (double *)malloc(size * sizeof(double)); }
	
	i = j = 0;
	while(!feof(file)){
		if(i >= size || j >= size) break;
		if(fscanf(file, "%lf", &val) != 1) return NULL;
		
		matrix[i][j] = val;
		
		j = (j + 1)%size;
		if(j == 0) i++;
	}
	
	fclose(file);	
	return matrix;
}

int main(){
	int N, i, j, e;
	char fname[50];
	double** matrix = NULL;

	matrix = read_matrix(fname, matrix, &N);
	
	if(!matrix){
		printf("The matrix couldn't be read!\n");
		exit(1);
	}
	
	return 0;
}

