#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include "util.h"

void ordenar_dados_com_MPI(int** v, int tamanho_v, int argc, char **argv);
void get_vetor_local(int *v, int tamanho_v, int **v_local, int *tamanho_v_local, int ranking, int num_proc);
void get_limites_vetor_local(int tamanho_v, int ranking, int num_proc, int *lim_inf, int *lim_sup);
void ordenar_dados(int **v, int tamanho_v, int ranking, int num_proc);
void mergeSort(int **v, int l, int r);
void merge(int **v, int l, int m, int r);
void fazer_merge_paralelo(int **v, int tamanho_v, int ranking, int num_proc);
void trocar_vetores_ordenados_localmente(int ranking, int vizinho, int* v, int tamanho_v,
                                         int **v2, int *tamanho_v2);
int* get_merge_vetores(int ranking, int dimensao, int *v, int tamanho_v, int* v2, int tamanho_v2);
int* obter_maiores_valores_vetores(int *v, int tamanho_v, int* v2, int tamanho_v2);
int* obter_menores_valores_vetores(int *v, int tamanho_v, int* v2, int tamanho_v2);
void get_vetor_ordenado(int **v, int tamanho_v, int *v_local, int tamanho_v_local, int ranking, int num_proc);

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
    int *v_local, tamanho_v_local;
    get_vetor_local(*v, tamanho_v, &v_local, &tamanho_v_local, ranking, num_proc);
    ordenar_dados(&v_local, tamanho_v_local, ranking, num_proc);


    fazer_merge_paralelo(&v_local, tamanho_v_local, ranking,num_proc);
    if(ranking==0)
        (*v)=v_local;

    if(ranking == 0){
        printf("tamanho_v: %d\nv: [", tamanho_v);
        for(int i=0; i < tamanho_v; i++){
            printf(" %d ", (*v)[i]);
        }
        printf("]\n");
    }

    MPI_Finalize();
    if(ranking != 0) exit(0);
}

void get_vetor_local(int *v, int tamanho_v, int **v_local, int *tamanho_v_local, int ranking, int num_proc){
    int lim_inf, lim_sup, *aux;
    get_limites_vetor_local(tamanho_v, ranking, num_proc, &lim_inf, &lim_sup);
    if(lim_inf == -1 || lim_sup == -1){
        *tamanho_v_local = 0;
        *v_local = NULL;
        return;
    }
    *tamanho_v_local = lim_sup - lim_inf + 1;
    aux = (int*) malloc(*tamanho_v_local*sizeof(int));
    for(int i=0, j=lim_inf; i < *tamanho_v_local; i++, j++){
        aux[i] = v[j];
    }
    *v_local = aux;
}

void get_limites_vetor_local(int tamanho_v, int ranking, int num_proc, int *lim_inf, int *lim_sup){
    double div = ceil((double)tamanho_v/ num_proc);
    *lim_inf = (div*ranking < tamanho_v) ? div*ranking : -1;
    *lim_sup = (div*(ranking+1)-1 < tamanho_v) ? div * (ranking+1) -1 :
                                                 (*lim_inf != -1) ? tamanho_v - 1 : -1;
}

void ordenar_dados(int **v, int tamanho_v, int ranking, int num_proc){
    //Código para fazer merge sort sequencial...
    mergeSort(v,0,tamanho_v-1);
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

void fazer_merge_paralelo(int **v, int tamanho_v, int ranking, int num_proc){
    int num_iteracoes = log(num_proc)/log(2);
    int vizinho, *v2, tamanho_v2, *aux;
    for(int i = 0; i < num_iteracoes; i++){
        if(ranking<num_proc){
            if (ranking>=num_proc/2)
            {
                vizinho = ranking-num_proc/2;
            }
            else{
            
                vizinho= num_proc/2+ranking;
            }
            trocar_vetores_ordenados_localmente(ranking, vizinho, *v, tamanho_v, &v2, &tamanho_v2);
            if (ranking<num_proc/2)
            {
                aux = get_merge_vetores(ranking, i, *v, tamanho_v, v2, tamanho_v2);
                free(*v);
                free(v2);
                *v = aux;
                tamanho_v+=tamanho_v2;
            }
            num_proc=num_proc/2;
        }
        else{
            return;
        }
        
    }
}

void trocar_vetores_ordenados_localmente(int ranking, int vizinho, int* v, int tamanho_v,
                                         int **v2, int *tamanho_v2){
    if(ranking > vizinho){
        MPI_Send(&tamanho_v, 1, MPI_INT, vizinho, 0, MPI_COMM_WORLD);
        MPI_Send(v, tamanho_v, MPI_INT, vizinho, 0, MPI_COMM_WORLD);
    }
    else{
        MPI_Recv(tamanho_v2, 1, MPI_INT, vizinho, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        *v2 = (int*) malloc(*tamanho_v2*sizeof(int));
        MPI_Recv(*v2, *tamanho_v2, MPI_INT, vizinho, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

int* get_merge_vetores(int ranking, int dimensao, int *v, int tamanho_v, int* v2, int tamanho_v2){
    int *aux;
    aux=obter_menores_valores_vetores(v,tamanho_v,v2,tamanho_v2);
    return aux;
}

int* obter_menores_valores_vetores(int *v, int tamanho_v, int* v2, int tamanho_v2){
    int *aux = (int*) malloc((tamanho_v+tamanho_v2)*sizeof(int));
    int j = 0, k = 0;
    for(int i = 0; i < tamanho_v+tamanho_v2; i++){
        if(j >= tamanho_v){
            aux[i] = v2[k];
            k++;
        }
        else if(k >= tamanho_v2){
            aux[i] = v[j];
            j++;
        }
        else if(v[j] <= v2[k]){
            aux[i] = v[j];
            j++;
        }
        else{
            aux[i] = v2[k];
            k++;
        }
    }
    return aux;
}