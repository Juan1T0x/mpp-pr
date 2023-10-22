#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <limits.h>

#define min 50
#define max 1000

double *generar_matriz_distancias(int n)
{
	double f;
	double *d = (double *) malloc(n * n *sizeof(double));

	srand(time(NULL) + getpid());
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			f = (double) rand() / ((double) RAND_MAX + 1);
			*(d + (i*n + j)) = (j != i) ? (min + f*(max - min)) : 0.0;
		}
	}
	return d;
}

void print_matrix(double *d, int n)
{
	for(int i = 0; i <= n*n; i++) {
		if ((i % n) != 0) { printf("%.0lf ", *(d + i)); }
		else { printf("\n%.0lf ", *(d + i)); }
	}
	printf("\n\n");
}

void print_solution(int n, const int *solucion, double valor)
{
	printf("\nSolution: ");
	for (int i = 0; i <= n; i++) {
		printf("%d ", solucion[i]);
	}
	printf("\nDistance: %.0lf\n", valor);
}

