#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include "util.h"

void ordenar_dados_com_MPI(int** v, int tamanho_v, int argc, char **argv);
void distribuir_dados();
void ordenar_dados_locais();
void get_limites_vetor(int tamanho_v, int ranking, int num_proc, int *lim_inf, int *lim_sup);
void mergeSort(int **v, int l, int r);
void merge(int **v, int l, int m, int r);
void fazer_merge_paralelo(int **v, int tamanho_v, int ranking, int num_proc);
void trocar_vetores_ordenados_localmente(int ranking, int vizinho, int* v, int tamanho_v,
                                         int **v2, int *tamanho_v2);
int* get_merge_vetores(int ranking, int dimensao, int *v, int tamanho_v, int* v2, int tamanho_v2);
int* obter_maiores_valores_vetores(int *v, int tamanho_v, int* v2, int tamanho_v2);
int* obter_menores_valores_vetores(int *v, int tamanho_v, int* v2, int tamanho_v2);


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
    /////////unir_dados_ordenados();
    MPI_Finalize();
    if(ranking != 0) exit(0);
}

void distribuir_dados(int **v,int tamanho_v,int ranking,int num_proc){
    //Código para distribuir os valores do vetor pelos processos...
    
}

void ordenar_dados_locais(int **v,int tamanho_v,int ranking,int num_proc){
    //Código para fazer merge sort sequencial...
    int l, r;
    get_limites_vetor(tamanho_v, ranking, num_proc, &l, &r);
    mergeSort(v,l,r);
}

void get_limites_vetor(int tamanho_v, int ranking, int num_proc, int *lim_inf, int *lim_sup){
    double div = ceil((double)tamanho_v/ num_proc);
    *lim_inf = div*ranking;
    *lim_sup = (ranking == num_proc - 1) ? tamanho_v-1 : div * (ranking+1) -1;
}

void mergeSort(int **v, int l, int r){
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

void merge(int **v, int l, int m, int r){
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

//**v deve ser o vetor inteiro para trabalharmos, não um pedaço dele
void fazer_merge_paralelo(int **v, int tamanho_v, int ranking, int num_proc){
    int num_iteracoes = log(num_proc)/log(2);
    int vizinho, *v2, tamanho_v2, *aux;
    for(int i = 0; i < num_iteracoes; i++){
        vizinho = ranking ^ (int) pow(2, i);
        trocar_vetores_ordenados_localmente(ranking, vizinho, *v, tamanho_v, &v2, &tamanho_v2);
        aux = get_merge_vetores(ranking, i, *v, tamanho_v, v2, tamanho_v2);
        free(*v);
        free(v2);
        *v = aux;
    }
}

void trocar_vetores_ordenados_localmente(int ranking, int vizinho, int* v, int tamanho_v,
                                         int **v2, int *tamanho_v2){
    if(ranking > vizinho){
        MPI_Send(&tamanho_v, 1, MPI_INT, vizinho, 0, MPI_COMM_WORLD);
        MPI_Send(v, tamanho_v, MPI_INT, vizinho, 0, MPI_COMM_WORLD);
        
        MPI_Recv(tamanho_v2, 1, MPI_INT, vizinho, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        *v2 = (int*) malloc(*tamanho_v2*sizeof(int));
        MPI_Recv(*v2, *tamanho_v2, MPI_INT, vizinho, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else{
        MPI_Recv(tamanho_v2, 1, MPI_INT, vizinho, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        *v2 = (int*) malloc(*tamanho_v2*sizeof(int));
        MPI_Recv(*v2, *tamanho_v2, MPI_INT, vizinho, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Send(&tamanho_v, 1, MPI_INT, vizinho, 0, MPI_COMM_WORLD);
        MPI_Send(v, tamanho_v, MPI_INT, vizinho, 0, MPI_COMM_WORLD);
    }
}

int* get_merge_vetores(int ranking, int dimensao, int *v, int tamanho_v, int* v2, int tamanho_v2){
    return (ranking & (int) pow(2, dimensao)) ? 
            obter_maiores_valores_vetores(v, tamanho_v, v2, tamanho_v2) : 
            obter_menores_valores_vetores(v, tamanho_v, v2, tamanho_v2);
}

int* obter_maiores_valores_vetores(int *v, int tamanho_v, int* v2, int tamanho_v2){
    int *aux = (int*) malloc(tamanho_v*sizeof(int));
    int j = tamanho_v2 - 1;
    for(int i = tamanho_v - 1; i <= 0; i--){
        if(v[i] >= v2[j]){
            aux[i] = v[i];
        }
        else{
            aux[i] = v2[j];
            j--;
        }
    }
    return aux;
}

int* obter_menores_valores_vetores(int *v, int tamanho_v, int* v2, int tamanho_v2){
    int *aux = (int*) malloc(tamanho_v*sizeof(int));
    int j = 0;
    for(int i = 0; i < tamanho_v; i++){
        if(v[i] <= v2[j]){
            aux[i] = v[i];
        }
        else{
            aux[i] = v2[j];
            j++;
        }
    }
    return aux;
}
