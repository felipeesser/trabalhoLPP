#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include "util.h"

void ordenar_dados_com_MPI(int** v, int tamanho_v, int argc, char **argv);
void get_vetor_local(int *v, int tamanho_v, int **v_local, int *tamanho_v_local, int ranking, int num_proc);
void get_limites_vetor_local(int tamanho_v, int ranking, int num_proc, int *lim_inf, int *lim_sup);
void ordenar_dados(int **v, int tamanho_v);
void mergeSort(int **v, int l, int r);
void merge(int **v, int l, int m, int r);
void fazer_merge_paralelo(int **v, int tamanho_v, int ranking, int num_proc);
void trocar_vetores_ordenados_localmente(int ranking, int vizinho, int* v, int tamanho_v,
                                         int **v2, int *tamanho_v2);
int* get_merge_vetores(int ranking, int dimensao, int *v, int tamanho_v, int* v2, int tamanho_v2);
void get_vetor_ordenado(int **v, int tamanho_v, int *v_local, int tamanho_v_local, int ranking, int num_proc);

int main(int argc, char **argv){
    char *arq_entrada = argv[1], *arq_saida = argv[2];
    int *v, tamanho_v;
    clock_t start, end;
    double cpu_time_used;
    
    //Obtendo o vetor a ser ordenado ("v")
    ler_entrada(arq_entrada, &v, &tamanho_v);
    
    start = clock();
    ordenar_dados_com_MPI(&v, tamanho_v, argc, argv);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    //Salvando o vetor ordenado ("v")
    escrever_saida(arq_saida, v, tamanho_v, cpu_time_used);

    liberar(&v);
    return 0;
}

void ordenar_dados_com_MPI(int** v, int tamanho_v, int argc, char **argv){
    //Iniciando MPI
    int ranking, num_proc;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ranking);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    int *v_local, tamanho_v_local;

    //Dividindo entre os processos o vetor original ("v"), em partes iguais.
    //O vetor que cada processo ordenará separadamente é "v_local"
    get_vetor_local(*v, tamanho_v, &v_local, &tamanho_v_local, ranking, num_proc);

    //Ordenando dados em "v_local"
    //Essa função realiza um merge sort sequencial
    ordenar_dados(&v_local, tamanho_v_local);

    //Unindo os vetores locais ("v_local") de cada processo para formar o vetor original ("v")
    //no processo mestre (ranking == 0)
    fazer_merge_paralelo(&v_local, tamanho_v_local, ranking,num_proc);
    if(ranking==0)
        (*v)=v_local;

    //Finalizando MPI e encerrado todos os processos exceto o que possui o vetor ordenado ("v")
    MPI_Finalize();
    
    if(ranking != 0) exit(0);
}

//Essa função recebe o vetor originalmente lido ("v") e retorna o vetor que o processo
//chamador deve ordenar ("v_local") sequencialmente
void get_vetor_local(int *v, int tamanho_v, int **v_local, int *tamanho_v_local, int ranking, int num_proc){
    //Obtendo os limites do vetor local do processo chamador ("v_local")
    int lim_inf, lim_sup, *aux;
    get_limites_vetor_local(tamanho_v, ranking, num_proc, &lim_inf, &lim_sup);
    if(lim_inf == -1 || lim_sup == -1){
        *tamanho_v_local = 0;
        *v_local = NULL;
        return;
    }

    //Alocando o vetor local do processo chamador ("v_local")
    *tamanho_v_local = lim_sup - lim_inf + 1;
    aux = (int*) malloc(*tamanho_v_local*sizeof(int));
    for(int i=0, j=lim_inf; i < *tamanho_v_local; i++, j++){
        aux[i] = v[j];
    }
    *v_local = aux;
}

//Essa função retorna os limites superior e inferior do vetor local ("v_local") do processo chamador
void get_limites_vetor_local(int tamanho_v, int ranking, int num_proc, int *lim_inf, int *lim_sup){
    double div = ceil((double)tamanho_v/ num_proc);
    *lim_inf = (div*ranking < tamanho_v) ? div*ranking : -1;
    *lim_sup = (div*(ranking+1)-1 < tamanho_v) ? div * (ranking+1) -1 :
                                                 (*lim_inf != -1) ? tamanho_v - 1 : -1;
}

//Essa função realiza um merge sort sequencial
void ordenar_dados(int **v, int tamanho_v){
    mergeSort(v,0,tamanho_v-1);
}

void mergeSort(int **v, int l, int r){
    if (l < r) {
        //(l+r)/2
        //Evita overflow para valores grandes de "l" e "r"
        int m = l + (r - l) / 2;  
        mergeSort(v, l, m);
        mergeSort(v, m + 1, r);
        merge(v, l, m, r);
    }
}

void merge(int **v, int l, int m, int r){
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
  
    //Vetores temporários
    int* L,*R;
    L=(int*)malloc(n1*sizeof(int));
    R=(int*)malloc(n2*sizeof(int));

  
    //Copiando porção do vetor original ("v") a ser ordenada para os vetores
    //temporários ("L" e "R")
    for (i = 0; i < n1; i++)
        L[i] = (*v)[l + i];
    for (j = 0; j < n2; j++)
        R[j] = (*v)[m + 1 + j];
  
    //Ordenando dados em "v"
    i = 0;
    j = 0;
    k = l;
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
  
    // Copiando dados de L e de R, caso algum sobre
    while (i < n1) {
        (*v)[k] = L[i];
        i++;
        k++;
    }
    while (j < n2) {
        (*v)[k] = R[j];
        j++;
        k++;
    }
    free(L);
    free(R);
}

//Essa função realiza o merge dos vetores locais ordenados em cada processo.
//Ao final de sua execução, o processo mestre (ranking == 0) possuirá o vetor original ("v")
//completamente ordenado
void fazer_merge_paralelo(int **v, int tamanho_v, int ranking, int num_proc){
    //A cada iteração do loop for, metade dos processos enviarão seus vetores locais ("v_local")
    //para a outra metade dos processos. Os processos que receberem os dados farão merge com seus
    //próprios vetores locais, formando um novo "v_local" com ambos os dados.
    //Apenas os processos que ainda não enviaram seus dados continuarão a execução, até que sobre
    //apenas o processo mestre (ranking == 0)
    int num_iteracoes = log(num_proc)/log(2);
    int vizinho, *v2, tamanho_v2, *aux;
    for(int i = 0; i < num_iteracoes; i++){
        if(ranking<num_proc){
            //Elegendo com que processo irá ocorrer a troca de dados
            if (ranking>=num_proc/2)
            {
                vizinho = ranking-num_proc/2;
            }
            else{
            
                vizinho= num_proc/2+ranking;
            }
            //Trocando dados entre os processos
            trocar_vetores_ordenados_localmente(ranking, vizinho, *v, tamanho_v, &v2, &tamanho_v2);
            if (ranking<num_proc/2)
            {
                //Fazendo merge dos dados recebidos
                aux = get_merge_vetores(ranking, i, *v, tamanho_v, v2, tamanho_v2);
                free(*v);
                free(v2);
                *v = aux;
                tamanho_v+=tamanho_v2;
            }
            num_proc=num_proc/2;
        }
        else{
            //eliminando processos que já enviaram seus dados
            MPI_Finalize();
            exit(0);
            return;
        }
        
    }
}

//Essa função realiza a troca de dados necessária para que o vetor local de um processo
//possa ser transmitido para outro processo
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

//Essa função realiza o merge de dois vetores ("v" e "v2") em um único vetor ("aux") e o retorna
int* get_merge_vetores(int ranking, int dimensao, int *v, int tamanho_v, int* v2, int tamanho_v2){
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
