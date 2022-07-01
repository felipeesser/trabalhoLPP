#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "util.h"
#include <omp.h>

void merge(int **V, int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
  
    int* L,*R;
    L=(int*)malloc(n1*sizeof(int));//elementos a esquerda da metade
    R=(int*)malloc(n2*sizeof(int));//elementos a direita da metade

    for (i = 0; i < n1; i++)
        L[i] = (*V)[l + i];
    for (j = 0; j < n2; j++)
        R[j] = (*V)[m + 1 + j];
  
    /* juncao ordenada das duas partes */
    i = 0; // indice de L[]
    j = 0; // indice de R[]
    k = l; // indice de V[] 
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            (*V)[k] = L[i];
            i++;
        }
        else {
            (*V)[k] = R[j];
            j++;
        }
        k++;
    }


    while (i < n1) {
        (*V)[k] = L[i];
        i++;
        k++;
    }
  
    while (j < n2) {
        (*V)[k] = R[j];
        j++;
        k++;
    }
    free(L);
    free(R);
}
  
void mergeSort(int **V, int l, int r)
{
    if (l < r) {
        int m = l + (r - l) / 2;
        #pragma omp parallel sections num_threads(2)
        {
            #pragma omp section
            {
mergeSort(V, l, m);
            }
            #pragma omp section
            {
mergeSort(V, m + 1, r);
            }
        }
        // divisao do vetor em duas partes
        
        //ordenacao 
        merge(V, l, m, r);
    }
}

int main(int argc, char **argv)
{
    int* V;//vetor a ser ordenado
    int N;//tmanho do vetor 
    double start,end,cpu_time_used;//contagem do tempo
    char *arq_entrada = argv[1];
    char *arq_saida = argv[2];
    ler_entrada(arq_entrada,&V,&N);
    start=omp_get_wtime();
    mergeSort(&V, 0,N - 1);
    end=omp_get_wtime();
    cpu_time_used=(end - start);
    escrever_saida(arq_saida,V,N,cpu_time_used);
    liberar(&V);
    
    return 0;
}