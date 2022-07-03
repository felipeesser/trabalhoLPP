#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "util.h"
#include <omp.h>
#include <math.h>

//Essa função realiza o merge de dois vetores em um vetor único. Os limites dos dois vetores são
//dados pelos parâmetros "l", "m" e "r".
void merge(int **V, int l, int m, int r)
{
    //Obtendo os tamanhos, "n1" e "n2", dos dois vetores que devemos fazer merge
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    //Criando os vetores auxiliares "L" e "R" e adicionando a eles os valores dos
    //dois vetores que devem ser mergedos
    int* L,*R;
    L=(int*)malloc(n1*sizeof(int));
    R=(int*)malloc(n2*sizeof(int));

    for (i = 0; i < n1; i++)
        L[i] = (*V)[l + i];
    for (j = 0; j < n2; j++)
        R[j] = (*V)[m + 1 + j];
  
    //Fazendo o merge ordenado dos dois vetores
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

    //Caso um dos vetores a serem mergeados, "L" ou "R", tenha acabado de ser adicionado
    //ao vetor resultante do merge, entretanto o outro vetor a ser mergeado ainda possuir valores
    //não adicionados ao vetor resultante, percorremos o vetor inacabado adicionando os
    //valores diretamente ao vetor resultante do merge
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

    //Apagando os vetores auxiliares "L" e "R"
    free(L);
    free(R);
}

//Essa função implementa o algoritmo Merge Sort sequencial
void mergeSort(int **V, int l, int r)
{
    if (l < r) {
        int m = l + (r - l) / 2;
        
        mergeSort(V, l, m);
                
        mergeSort(V, m + 1, r);
       
        merge(V, l, m, r);
    }
}

//Essa função realiza toda a ordenação do vetor pelo algoritmo Merge Sort paralelo
void paralel_merge(int **V, int size,int nthreads){
    int div;
    int l,r,m;
    int iterate=nthreads;
    int aux=1;
    #pragma omp parallel num_threads(nthreads) private(l,r,m)
    {
        //Dividindo o vetor original em partes iguais para que cada thread ordene uma delas
        //utilizando o algoritmo Merge Sort sequencial.
        //Para obtermos a partição que cabe a cada thread ordenar, determinamos os limites
        //da mesma (variáveis "l" e "r").
        div=size/omp_get_num_threads();
        l=(div* (omp_get_thread_num()));
        if(omp_get_thread_num()==(omp_get_num_threads()-1))
            r=size-1;
        else
            r=( div * (omp_get_thread_num()+1))-1;
        
        //Ordenando sequencialmente uma partição do vetor
        mergeSort(V,l,r);

        //Sincronizando todas as threads a fim de ter certeza de que, antes de iniciarmos o merge,
        //Todas as threads já ordenaram sequencialmente suas partições do vetor original
        #pragma omp barrier

        //Realizando o merge das partições do vetor ordenadas separadamente
        //O merge entre as partições é feito a cada iteração, sempre entre pares de partições consecutivas
        //Ex.: Supondo 8 partições, a execução seria a seguinte:
        //(início) T0, T1, T2, T3, T4, T5, T6, T7 -> (1ª iteração) T0T1, T2T3, T4T5, T6T7 ->
        // -> (2ª iteração) T0T1T2T3, T4T5T6T7 -> (3ª iteração) T0T1T2T3T4T5T6T7
        //Para isso, utilizamos cálculos baseados em potências de 2 de acordo com o número da 
        //iteração. Ex.: 1ª iteração -> 2^1, 2ª iteração -> 2^2, 3ª iteração -> 2^3, etc.
        while (iterate!=1)
        {
            //A thread que realizará o merge deve ser múltipla de 2^(iteração atual)
            if(omp_get_thread_num()%((int)pow(2,aux))==0){
                m=r;
                if(omp_get_thread_num()+(int)pow(2,aux)==omp_get_num_threads())
                    r=size-1;
                else
                    r=( div * (omp_get_thread_num()+(int)pow(2,aux)))-1;
                merge(V,l,m,r);
            }

            //Necessitamos utilizar as duas barreiras abaixo porque a condição de parada do loop
            //"while" depende do valor da variável "iterate" que é compartilhada entre todas as
            //threads e necessita ser atualizada a cada iteração do programa.
            //Assim, a primeira barreira garante que todas as threads executaram a iteração atual
            //do ciclo "while", ou seja, elas passaram pela condição de parada "while"
            //A segunda barreira garante que a variável "iterate" foi atualizada pela thread 0
            //antes que as demais threads verificassem se executariam por mais um loop ou não.
            #pragma omp barrier
            //Atualizando as variáveis compartilhadas para a próxima iteração
            if (omp_get_thread_num()==0){
                iterate=iterate/2;
                aux+=1;
            }
            #pragma omp barrier
        } 
    }
    
}

int main(int argc, char **argv)
{
    int* V; //vetor a ser ordenado
    int N; //tmanho do vetor 
    double start,end,cpu_time_used; //contagem do tempo
    char *arq_entrada = argv[1];
    char *arq_saida = argv[2];
    int nthreads = atoi(argv[3]);
    ler_entrada(arq_entrada,&V,&N);
    start=omp_get_wtime();
    paralel_merge(&V,N,nthreads);
    end=omp_get_wtime();
    cpu_time_used=(end - start);
    escrever_saida(arq_saida,V,N,cpu_time_used);
    liberar(&V);
    return 0;
}