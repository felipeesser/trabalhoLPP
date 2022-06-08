#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include "util.h"

void ordenar_dados_com_MPI(int** v, int tamanho_v, int argc, char **argv);
void distribuir_dados();
void ordenar_dados_locais();
void unir_dados_ordenados();

int main(int argc, char **argv){
    char *arq_entrada = argv[1], *arq_saida = argv[2];
    int *v, tamanho_v;
    clock_t start, end;
    double cpu_time_used;
    ler_entrada(arq_entrada, &v, &tamanho_v);
    start = clock();
    ordenar_dados_com_MPI(&v, tamanho_v, argc, argv);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    escrever_saida(arq_saida, v, tamanho_v, cpu_time_used);
    liberar(&v);
    return 0;
}

void ordenar_dados_com_MPI(int** v, int tamanho_v, int argc, char **argv){
    int ranking, num_proc;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ranking);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    distribuir_dados();
    ordenar_dados_locais(v,tamanho_v,ranking,num_proc); 
    unir_dados_ordenados();
    MPI_Finalize();
    if(ranking != 0) exit(0);
}

void merge(int **v, int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
  
    /* create temp arrays */
    int* L,*R;
    L=(int*)malloc(n1*sizeof(int));
    R=(int*)malloc(n2*sizeof(int));

  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = (*v)[l + i];
    for (j = 0; j < n2; j++)
        R[j] = (*v)[m + 1 + j];
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            (*v)[k] = L[i];
            i++;
        }
        else {
            (*v)[k] = R[j];
            j++;
        }
        k++;
    }
  
    /* Copy the remaining elements of L[], if there
    are any */
    while (i < n1) {
        (*v)[k] = L[i];
        i++;
        k++;
    }
  
    /* Copy the remaining elements of R[], if there
    are any */
    while (j < n2) {
        (*v)[k] = R[j];
        j++;
        k++;
    }
    free(L);
    free(R);
}

void mergeSort(int **v, int l, int r)
{
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l + (r - l) / 2;
  
        // Sort first and second halves
        mergeSort(v, l, m);
        mergeSort(v, m + 1, r);
  
        merge(v, l, m, r);

    }
}


void ordenar_dados_locais(int **v,int tamanho_v,int ranking,int num_proc){
    //Código para fazer merge sort sequencial...
    double div=ceil((double)tamanho_v/ num_proc);
    int l=(div* (ranking));
    int r;
    if(ranking==(num_proc-1))
        r=tamanho_v-1;
    else
        r=( div * (ranking+1))-1;
    mergeSort(v,l,r);
}

void distribuir_dados(int **v,int tamanho_v,int ranking,int num_proc){
    //Código para distribuir os valores do vetor pelos processos...
    
}

void unir_dados_ordenados(){
    //Código para fazer comunicação via hipercubo...
}
