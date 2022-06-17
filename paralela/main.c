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
void fazer_merge_paralelo(int **v, int tamanho_v, int ranking, int num_proc, int num_proc_ordenados);
int get_vizinho(int ranking, int fator, int num_proc_ordenados, int processo_principal, int prim_proc, int ultimo_proc);
void trocar_vetores_ordenados_localmente(int ranking, int vizinho, int* v, int tamanho_v,
                                         int **v2, int *tamanho_v2);
int* get_merge_vetores(int ranking, int vizinho, int *v, int tamanho_v, int* v2, int tamanho_v2);
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


    int num_iteracoes = log(num_proc)/log(2);
    for(int i=0; i < num_iteracoes; i++){
        fazer_merge_paralelo(&v_local, tamanho_v_local, ranking, num_proc, pow(2, i));


        /***********
        printf("(rank %d) tamanho_v_local: %d\nv_local: [", ranking, tamanho_v_local);
        for(int i=0; i < tamanho_v_local; i++){
            printf("(rank %d) %d ", ranking, v_local[i]);
        }
        printf("]\n");
        *****************///////////////
    }
    

    get_vetor_ordenado(v, tamanho_v, v_local, tamanho_v_local, ranking, num_proc);

    
    if(ranking == 0){
        printf("tamanho_v: %d\nv: [", tamanho_v);
        for(int i=0; i < tamanho_v; i++){
            printf(" %d ", (*v)[i]);
        }
        printf("]\n");
    }
    

    if(v_local != NULL) free(v_local);
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

//dim_proc_ordenados = log2(número de processos ordenados)
//Ex.: Se há 8 processos ordenados, então dim_proc_ordenados = 3
void fazer_merge_paralelo(int **v, int tamanho_v, int ranking, int num_proc, int num_proc_ordenados){
    int dim_proc_ordenados = log(num_proc_ordenados)/log(2);
    //O prim_proc é o processo com menor ranking que está em um dos 2 grupos que serão combinados
    //O ultimo_proc é o processo com maior ranking que está em um dos 2 grupos que serão combinados
    int prim_proc = ranking - (ranking % (2*num_proc_ordenados));
    int ultimo_proc = prim_proc + 2*num_proc_ordenados - 1;
    int processo_principal = prim_proc + dim_proc_ordenados;
    int vizinho, *v2, tamanho_v2, *aux;

    int fator = 0;
    for(int i = dim_proc_ordenados, j = 0; i >= 0; i--, j++){
        fator = pow(2, i);
        if(j % 2 == 0) fator = -fator;
        vizinho = get_vizinho(ranking, fator, num_proc_ordenados, processo_principal, prim_proc, ultimo_proc);

        

        printf("(rank %d) fator %d vizinho %d AUX: [\n", ranking, fator, vizinho);





        trocar_vetores_ordenados_localmente(ranking, vizinho, *v, tamanho_v, &v2, &tamanho_v2);
        aux = get_merge_vetores(ranking, vizinho, *v, tamanho_v, v2, tamanho_v2);

        
        /*
        for(int i=0; i < tamanho_v; i++){
            printf("(rank %d AUX) %d ", ranking, aux[i]);
        }
        printf("]\n");
        */

        
        free(*v);
        free(v2);
        *v = aux;
    }
}

int get_vizinho(int ranking, int fator, int num_proc_ordenados, int processo_principal, int prim_proc, int ultimo_proc){
    int proc_escolhidos[2*num_proc_ordenados], num_proc_escolhidos=0;
    for(int i=0; i < 2*num_proc_ordenados; i++)
        proc_escolhidos[i] = -1;

    int proc_atual = processo_principal, vizinho_atual;
    int proc_atual_invalido = 0;
    while(1){
        vizinho_atual = proc_atual + fator;
        if(vizinho_atual < prim_proc){
            vizinho_atual = ultimo_proc - (prim_proc - vizinho_atual) + 1;
        }
        else if(vizinho_atual > ultimo_proc){
            vizinho_atual = prim_proc + (vizinho_atual - ultimo_proc) - 1;
        }

        for(int i = 0; i < num_proc_escolhidos; i++){
            if(vizinho_atual == proc_escolhidos[i]){
                proc_atual_invalido = 1;
                break;
            }
        }

        if(!proc_atual_invalido){
            if(proc_atual == ranking) return vizinho_atual;
            if(vizinho_atual == ranking) return proc_atual;

            proc_escolhidos[num_proc_escolhidos] = proc_atual;
            proc_escolhidos[num_proc_escolhidos+1] = vizinho_atual;
            num_proc_escolhidos += 2;

            proc_atual_invalido = 1;
        }

        while(proc_atual_invalido){
            proc_atual += 1;
            if(proc_atual > ultimo_proc){
                proc_atual = prim_proc;
            }
            proc_atual_invalido = 0;
            for(int i = 0; i < num_proc_escolhidos; i++){
                if(proc_atual == proc_escolhidos[i]){
                    proc_atual_invalido = 1;
                    break;
                }
            }
        }
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

int* get_merge_vetores(int ranking, int vizinho, int *v, int tamanho_v, int* v2, int tamanho_v2){
    return (ranking > vizinho) ? 
            obter_maiores_valores_vetores(v, tamanho_v, v2, tamanho_v2) : 
            obter_menores_valores_vetores(v, tamanho_v, v2, tamanho_v2);
}

int* obter_maiores_valores_vetores(int *v, int tamanho_v, int* v2, int tamanho_v2){
    int *aux = (int*) malloc(tamanho_v*sizeof(int));
    int j = tamanho_v - 1, k = tamanho_v2 - 1;
    for(int i = tamanho_v - 1; i >= 0; i--){
        if(j < 0){
            aux[i] = v2[k];
            k--;
        }
        else if(k < 0){
            aux[i] = v[j];
            j--;
        }
        else if(v[j] >= v2[k]){
            aux[i] = v[j];
            j--;
        }
        else{
            aux[i] = v2[k];
            k--;
        }
    }
    return aux;
}

int* obter_menores_valores_vetores(int *v, int tamanho_v, int* v2, int tamanho_v2){
    int *aux = (int*) malloc(tamanho_v*sizeof(int));
    int j = 0, k = 0;
    for(int i = 0; i < tamanho_v; i++){
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

void get_vetor_ordenado(int **v, int tamanho_v, int *v_local, int tamanho_v_local, int ranking, int num_proc){
    MPI_Gather(v_local, tamanho_v_local, MPI_INT, *v, tamanho_v_local, MPI_INT, 0, MPI_COMM_WORLD);
}