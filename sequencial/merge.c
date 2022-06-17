#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "util.h"
#include <string.h>

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
  
        // divisao do vetor em duas partes
        mergeSort(V, l, m);
        mergeSort(V, m + 1, r);

        //ordenacao 
        merge(V, l, m, r);
    }
}

int main(int argc, char **argv)
{
    FILE* f;
    int* V;//vetor a ser ordenado
    int N;//tmanho do vetor 
    double start,end,cpu_time_used;//contagem do tempo
    int num_trials=atoi(argv[1]);
    char *arq_saida = argv[2];
    char nome[30]="../aleatorio/entradas/entrada";
    char num_trial[50];
    char aux[80];

    for (int i = 1; i <= num_trials; i++)
    {   
        strcpy(aux,nome);
        sprintf(num_trial,"%d",i);
        ler_entrada(strcat(aux,num_trial),&V,&N);
        start=clock();
        mergeSort(&V, 0,N - 1);
        end=clock();
        cpu_time_used=((double) (end - start)) / CLOCKS_PER_SEC;
        escrever_saida(arq_saida,V,N,cpu_time_used);
        liberar(&V);
    }
    
    return 0;
}