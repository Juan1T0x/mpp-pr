#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "../include/ga.h"

#define PRINT 1

unsigned int seed = 1;


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

double aplicar_ga(const double *d, int n, int n_gen, int tam_pob, int *sol, double m_rate)
{
	int i, g, mutation_start;

	/*
	double fitness_anterior = -1.0;
	double threshold = 0.1;
	double improvement = 100.0;
	*/
	
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
			mutar(poblacion[i], n, m_rate);
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
		/*
		// Verificar si el fitness actual es al menos un N% mejor que el anterior
		improvement = ((fitness_anterior - poblacion[0]->fitness) / fitness_anterior) * 100.0;

        if (improvement <= threshold)
        {
            printf("Condición de parada cumplida, el valor del fitness actual no es mejor que el anterior según el criterio indicado\n");
            break;
        }

        fitness_anterior = poblacion[0]->fitness;
		*/

	}
	
	memmove(sol, poblacion[0]->array_int, n*sizeof(int));
	
	// almacena el mejor valor obtenido para el fitness
	double value = (poblacion[0]->fitness);
	
	// Liberar memoria de cada individuo y su array_int
	for (i = 0; i < tam_pob; i++) {
	    free(poblacion[i]->array_int);
	    free(poblacion[i]);
	}

	// Liberar memoria del array de punteros a Individuo
	free(poblacion);
	
	// devuelve el valor obtenido para el fitness
	return value;
}

void cruzar(Individuo *padre1, Individuo *padre2, Individuo *hijo1, Individuo *hijo2, int n)
{
	// Elegir un punto (o puntos) de corte aleatorio a partir del que se realiza el intercambio de los genes. 
    int punto = aleatorio(n-1) + 1; // Selecciona punto de corte (evitando el 0)
    
    // Copia genes hasta el punto de corte
    for(int i=1; i < punto; i++) {
        hijo1->array_int[i] = padre1->array_int[i];
        hijo2->array_int[i] = padre2->array_int[i];
    }

    // Copia genes después del punto de corte evitando duplicados
    for(int i=punto; i<n; i++) {
        if(!search_element(hijo1->array_int, i, padre2->array_int[i])) {
            hijo1->array_int[i] = padre2->array_int[i];
        } else {
            for(int j=1; j<n; j++) {
                if(!search_element(hijo1->array_int, i, padre2->array_int[j])) {
                    hijo1->array_int[i] = padre2->array_int[j];
                    break;
                }
            }
        }
        
        if(!search_element(hijo2->array_int, i, padre1->array_int[i])) {
            hijo2->array_int[i] = padre1->array_int[i];
        } else {
            for(int j=1; j<n; j++) {
                if(!search_element(hijo2->array_int, i, padre1->array_int[j])) {
                    hijo2->array_int[i] = padre1->array_int[j];
                    break;
                }
            }
        }
    }
}

void invertir(int *a, int i, int j)
{
    while (i < j) {
        int temp = a[i];
        a[i] = a[j];
        a[j] = temp;
        i++;
        j--;
    }
}

void mutar(Individuo *actual, int n, double m_rate)
{
    //double m_rate = 0.15; // Modificar valor si conviene

    for (int m = 0; m < n * m_rate; m++) {
        int i = aleatorio(n-1) + 1; // Asegura que i > 0
        int j = aleatorio(n-1) + 1;

        if (i > j) {
            int temp = i;
            i = j;
            j = temp;
        }

        invertir(actual->array_int, i, j);
    }
}

double distancia_ij(const double *d, int i, int j, int n)
{
    return d[i * n + j];
}

void fitness(const double *d, Individuo *individuo, int n)
{
    double fit = 0.0;
    for(int i=0; i<n-1; i++) {
        fit += distancia_ij(d, individuo->array_int[i], individuo->array_int[i+1], n);
    }
    individuo->fitness = fit;
}
