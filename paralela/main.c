#include <stdio.h>
#include <time.h>
#include <mpi.h>

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
    ordenar_dados_com_MPI(v, tamanho_v, argc, argv);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    escrever_saida(arq_saida, v, cpu_time_used);
    
    return 0;
}

void ordenar_dados_com_MPI(int** v, int tamanho_v, int argc, char **argv){
    int ranking, num_proc;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ranking);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    if(ranking == 0) distribuir_dados();
    ordenar_dados_locais();
    unir_dados_ordenados();
    MPI_Finalize();
}

void distribuir_dados(){
    //Código para distribuir os valores do vetor pelos processos...
}

void ordenar_dados_locais(){
    //Código para fazer merge sort sequencial...
}

void unir_dados_ordenados(){
    //Código para fazer comunicação via hipercubo...
}
