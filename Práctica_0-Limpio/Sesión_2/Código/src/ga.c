#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "../include/ga.h"

#define PRINT 1

int aleatorio(int n) {
	return (rand() % n);  // genera un numero aleatorio entre 0 y n-1
}

int search_element(int *array, int end, int element)
{
	int i=0;
	int found=0;
	
	// comprueba que un elemento no está incluido en el individuo (el cual no admite enteros repetidos)
	while((i < end) && ! found) {
		if(array[i] == element) {
			found = 1;
		}
		i++;
	}
        return found;
}

int find_element(int *array, int end, int element)
{
        int pos = 0;
	for(int i = 0; i < end; i++) {
             if(array[i] == element) {
                 pos = i;
                 break;
             }
        }
        return pos; // Posición del elemento encontrado
}

int *crear_individuo(int n)
{
        // El primer elemento del individuo siempre será el 0, por ejemplo.
	int i=1, value;
	int *individuo = (int *) malloc(n * sizeof(int));
	
	// inicializa array de elementos
	memset(individuo, 0, n * sizeof(int));
	
	while(i < n) {
		value = aleatorio(n);
		// si el nuevo elemento no está en el array...
		if(!search_element(individuo, i, value)) {
			individuo[i] = value;  // lo incluimos
			i++;
		}
	}
	return individuo;
}

int comp_fitness(const void *a, const void *b) {
	/* qsort pasa un puntero al elemento que está ordenando */
	return (*(Individuo **)a)->fitness - (*(Individuo **)b)->fitness;
}

double aplicar_ga(const double *d, int n, int n_gen, int tam_pob, int *sol)
{
	int i, g, mutation_start;
	
	// crea poblacion inicial (array de individuos)
	Individuo **poblacion = (Individuo **) malloc(tam_pob * sizeof(Individuo *));
	assert(poblacion);
	
	// crea cada individuo (array de enteros aleatorios)
	for(i = 0; i < tam_pob; i++) {
		poblacion[i] = (Individuo *) malloc(sizeof(Individuo));
		poblacion[i]->array_int = crear_individuo(n);
		
		// calcula el fitness del individuo
		fitness(d, poblacion[i], n);
	}
	
	// ordena individuos segun la funcion de bondad (menor "fitness" --> mas aptos)
	qsort(poblacion, tam_pob, sizeof(Individuo *), comp_fitness);
	
	// evoluciona la poblacion durante un numero de generaciones
	for(g = 0; g < n_gen; g++)
	{
		// los hijos de los ascendientes mas aptos sustituyen a la ultima mitad de los individuos menos aptos
		for(i = 0; i < (tam_pob/2) - 1; i += 2) {
			cruzar(poblacion[i], poblacion[i+1], poblacion[tam_pob/2 + i], poblacion[tam_pob/2 + i + 1], n);
		}
		
		// por ejemplo, inicia la mutacion a partir de 1/4 de la poblacion.
                // puede haber otras opciones pero dejar el primer individuo sin modificar siempre
		mutation_start = tam_pob/4;
		
		// muta 3/4 partes de la poblacion
		for(i = mutation_start; i < tam_pob; i++) {
			mutar(poblacion[i], n);
		}
		
		// recalcula el fitness del individuo
		for(i = 0; i < tam_pob; i++) {
			fitness(d, poblacion[i], n);
		}
		
		// ordena individuos segun la funcion de bondad (menor "fitness" --> mas aptos)
		qsort(poblacion, tam_pob, sizeof(Individuo *), comp_fitness);
		
		if (PRINT) {
			printf("Generacion %d - ", g);
			printf("Fitness = %.0lf\n", (poblacion[0]->fitness));
		}
	}
	
	memmove(sol, poblacion[0]->array_int, n*sizeof(int));
	
	// almacena el mejor valor obtenido para el fitness
	double value = (poblacion[0]->fitness);
	
	// se libera la memoria reservada
	free(poblacion);
	
	// devuelve el valor obtenido para el fitness
	return value;
}

void cruzar(Individuo *padre1, Individuo *padre2, Individuo *hijo1, Individuo *hijo2, int n)
{
	// Elegir un punto (o puntos) de corte aleatorio a partir del que se realiza el intercambio de los genes. 
	
	// Entonces, por ejemplo, los primeros genes del padre1 van al hijo1, y los primeros del padre2 al hijo2.
        // Se debe evitar en cada paso la introduccion de duplicados en los hijos
	// Y los restantes genes de cada hijo son del otro padre, respectivamente.

        // Otra opción podría ser factibilizar a posteriori, despues de generar los descendientes: eliminar posibles 
        // repetidos de ambos hijos. Si encuentro algún elemento repetido en el hijo, lo cambio por otro que no este el array
}

void invertir(int *a, int k)
{
        int t;
	// Uno por uno invierte los elementos de a[0..k-1]
}

void mutar(Individuo *actual, int n)
{
        double m_rate = 0.15; // modificar valor si conviene

	// Implementación recomendada (aunque puede ser cualquier otra que se considere adecuada para este problema): 
	// Reverse Sequence Mutation (RSM), donde elegimos una secuencia S limitada por dos posiciones i, j
        // elegidas aleatoriamente con i<j, e i>0 para no modificar nunca el 1er elemento. El orden de los elementos en 
	// esta secuencia será invertido, por ejemplo con i=1, j=4: (1,2,3,4,5,6) --> (1,5,4,3,2,6).

        // Usar la variable m_rate para establecer la intensidad (iteraciones) de la mutación, teniendo en cuenta que
	// si el valor es demasiado pequeño la convergencia es muy pequeña y si es demasiado puede diverger.
}

double distancia_ij(const double *d, int i, int j, int n)
{
	double dist = 0.0;
	
	// Devuelve la distancia entre dos elementos de la matriz 'd'
	return dist;
}

void fitness(const double *d, Individuo *individuo, int n)
{
	// Determina la calidad del individuo calculando la suma de la distancia entre cada par de ciudades consecutivas en el array
}
