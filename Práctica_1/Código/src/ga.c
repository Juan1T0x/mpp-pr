#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <omp.h>
#include <time.h>

#include "../include/ga.h"

#define PRINT 1

#define N_THREADS 4
#define CHUNKSIZE 2

unsigned int seed = 1;
#pragma omp threadprivate(seed)

int aleatorio(int n) {
	//return (rand() % n);  // genera un numero aleatorio entre 0 y n-1
	return (rand_r(&seed) % n);  // genera un numero aleatorio entre 0 y n-1 (thread-safe)
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

	/*
	double fitness_anterior = -1.0;
	double threshold = 0.1;
	double improvement = 100.0;
	*/
	
	// crea poblacion inicial (array de individuos)
	Individuo **poblacion = (Individuo **) malloc(tam_pob * sizeof(Individuo *));
	assert(poblacion);

	#pragma omp parallel
	{
		seed = time(NULL) ^ omp_get_thread_num();
		#pragma omp for
		for(i = 0; i < tam_pob; i++) {
			poblacion[i] = (Individuo *) malloc(sizeof(Individuo));
			poblacion[i]->array_int = crear_individuo(n);
			
			// calcula el fitness del individuo
			fitness(d, poblacion[i], n);
		}
	}
	
	/*
	// crea cada individuo (array de enteros aleatorios)
	for(i = 0; i < tam_pob; i++) {
		poblacion[i] = (Individuo *) malloc(sizeof(Individuo));
		poblacion[i]->array_int = crear_individuo(n);
		
		// calcula el fitness del individuo
		fitness(d, poblacion[i], n);
	}
	*/
	
	// ordena individuos segun la funcion de bondad (menor "fitness" --> mas aptos)
	qsort(poblacion, tam_pob, sizeof(Individuo *), comp_fitness);
	
	// evoluciona la poblacion durante un numero de generaciones
	for(g = 0; g < n_gen; g++)
	{
		// los hijos de los ascendientes mas aptos sustituyen a la ultima mitad de los individuos menos aptos
		for(i = 0; i < (tam_pob/2) - 1; i += 2) {
			cruzar_paralel(poblacion[i], poblacion[i+1], poblacion[tam_pob/2 + i], poblacion[tam_pob/2 + i + 1], n);
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
	int i;
	//#pragma omp parallel for private(i) schedule(static, CHUNKSIZE) num_threads(N_THREADS)
	//#pragma omp parallel for private(i) schedule(dynamic, CHUNKSIZE) num_threads(N_THREADS)
	//#pragma omp parallel for private(i) schedule(guided, CHUNKSIZE) num_threads(N_THREADS)
    for(i=1; i < punto; i++) {
        hijo1->array_int[i] = padre1->array_int[i];
        hijo2->array_int[i] = padre2->array_int[i];
    }

    // Copia genes después del punto de corte evitando duplicados
	//#pragma omp parallel for private(i) schedule(static, CHUNKSIZE) num_threads(N_THREADS)
	//#pragma omp parallel for private(i) schedule(dynamic, CHUNKSIZE) num_threads(N_THREADS)
	//#pragma omp parallel for private(i) schedule(guided, CHUNKSIZE) num_threads(N_THREADS)
    for(i=punto; i<n; i++) {
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

void cruzar_paralel(Individuo *padre1, Individuo *padre2, Individuo *hijo1, Individuo *hijo2, int n){
	
	omp_set_num_threads(N_THREADS);

	int punto = aleatorio(n-1) + 1; // Selecciona punto de corte (evitando el 0)

    #pragma omp parallel
    {
        #pragma omp sections nowait
        {
            #pragma omp section
            {
                for(int i=1; i < punto; i++) {
                    hijo1->array_int[i] = padre1->array_int[i];
                    hijo2->array_int[i] = padre2->array_int[i];
                }
            }

            #pragma omp section
            {
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

void mutar(Individuo *actual, int n)
{
    double m_rate = 0.15; // Modificar valor si conviene
    int num_iterations = (int)(n * m_rate); // Calcular el número de iteraciones como entero

    int m;
    //#pragma omp parallel for private(m) schedule(static, CHUNKSIZE) num_threads(N_THREADS)
	//#pragma omp parallel for private(m) schedule(dynamic, CHUNKSIZE) num_threads(N_THREADS)
	//#pragma omp parallel for private(m) schedule(guided, CHUNKSIZE) num_threads(N_THREADS)
    for (m = 0; m < num_iterations; m++) {
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
	int i;
	//#pragma omp parallel for private(i) schedule(static, CHUNKSIZE) num_threads(N_THREADS)
	//#pragma omp parallel for private(i) schedule(dynamic, CHUNKSIZE) num_threads(N_THREADS) reduction(+:fit)
	//#pragma omp parallel for private(i) schedule(guided, CHUNKSIZE) num_threads(N_THREADS) reduction(+:fit)
    for(i=0; i<n-1; i++) {
        fit += distancia_ij(d, individuo->array_int[i], individuo->array_int[i+1], n);
    }
    individuo->fitness = fit;
}

void fitness_critical(const double *d, Individuo *individuo, int n)
{
	omp_set_num_threads(N_THREADS);

	double fit = 0.0;
	double temp;
	individuo->fitness = 0.0;
	#pragma omp parallel for private(temp) shared(fit)
	for(int i=0; i<n-1; i++) {
		temp = distancia_ij(d, individuo->array_int[i], individuo->array_int[i+1], n);
		#pragma omp critical
		fit += temp;
	}
	individuo->fitness = fit;
}

void fitness_atomic(const double *d, Individuo *individuo, int n)
{
	omp_set_num_threads(N_THREADS);

	double fit = 0.0;
	double temp;
	individuo->fitness = 0.0;
	#pragma omp parallel for private(temp) shared(fit)
	for(int i=0; i<n-1; i++) {
		temp = distancia_ij(d, individuo->array_int[i], individuo->array_int[i+1], n);
		#pragma omp atomic
		fit += temp;
	}
	individuo->fitness = fit;
}

void fitness_reduction(const double *d, Individuo *individuo, int n)
{

	omp_set_num_threads(N_THREADS);

	double fit = 0.0;
	double temp;
	individuo->fitness = 0.0;

	#pragma omp parallel for reduction(+:fit) private(temp)
		for(int i=0; i<n-1; i++) {
		temp = distancia_ij(d, individuo->array_int[i], individuo->array_int[i+1], n);
		fit += temp;
	}
	individuo->fitness = fit;
}

void mezclar(Individuo **poblacion, int izq, int med, int der)
{	
	int i, j, k;
	
	Individuo **pob = (Individuo **) malloc((der - izq)*sizeof(Individuo *));
	assert(pob);

	for(i = 0; i < (der - izq); i++) {
		pob[i] = (Individuo *) malloc(sizeof(Individuo));
	}
	
	k = 0;
	i = izq;
	j = med;
	while( (i < med) && (j < der) ) {
		if (poblacion[i]->fitness > poblacion[j]->fitness) {
			// copiar poblacion[i++] en pob[k++]
		}
		else {
			// copiar poblacion[j++] en pob[k++]
		}
	}
	
	for(; i < med; i++) {
		// copiar poblacion[i] en pob[k++]
	}
	
	for(; j < der; j++) {
		// copiar poblacion[j] en pob[k++]
	}
	
	i = 0;
	while(i < (der - izq)) {
		// copiar pob[i] en poblacion[i + izq]
		// ...
		// liberar individuo 'i'
		free(pob[i]);
		i++;
	}
	free(pob);
}

void mergeSort(Individuo **poblacion, int izq, int der) 
{ 
	int med = (izq + der)/2; 
	if ((der - izq) < 2) { return; } 
  	mergeSort(poblacion, izq, med); 
	mergeSort(poblacion, med, der);   
	mezclar(poblacion, izq, med, der); 
}
