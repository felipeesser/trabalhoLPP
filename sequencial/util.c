#include "util.h"

void ler_entrada(char *nome_arq, int **v, int *tamanho_v){
    FILE* arq = fopen(nome_arq, "rt");
    fscanf(arq, "%d", tamanho_v);
    *v= (int*) malloc((*tamanho_v)*sizeof(int));
    for(int i = 0; i < *tamanho_v; i++) fscanf(arq, "%d", &(*v)[i]);
    fclose(arq);
}

void escrever_saida(char *nome_arq, int *v, int tamanho_v, double tempo){
    FILE* arq = fopen(nome_arq, "at");
    fprintf(arq, "Vetor Ordenado: [");
    for(int i = 0; i < tamanho_v; i++){
        fprintf(arq, "%d", v[i]);
        if(i != tamanho_v - 1) fprintf(arq, ", ");
    }
    fprintf(arq, "]\n");
    fprintf(arq, "Tempo de CPU: %lf\n", tempo);
    fclose(arq);
    FILE* arq_temp = fopen("tempo_seq", "w");
    fprintf(arq_temp, "%d %lf\n", tamanho_v,tempo);
    fclose(arq_temp);
}

void liberar(int **v){
    free(*v);
}
