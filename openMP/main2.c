#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "util.h"
#include <omp.h>
#include <math.h>

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
        
        mergeSort(V, l, m);
                
        mergeSort(V, m + 1, r);
       
        merge(V, l, m, r);
    }
}

void paralel_merge(int **V, int size,int nthreads){
    int div;
    int l,r,m;
    int iterate=nthreads;
    int aux=1;
    #pragma omp parallel num_threads(nthreads) private(l,r,m)
    {
        //ordenacao dos vetores locais 
        div=size/omp_get_num_threads();
        l=(div* (omp_get_thread_num()));
        if(omp_get_thread_num()==(omp_get_num_threads()-1))
            r=size-1;
        else
            r=( div * (omp_get_thread_num()+1))-1;
        mergeSort(V,l,r);
        #pragma omp barrier
        //fim ordenacao vetores locais
        //merge dos vetores locais
        while (iterate!=1)
        {
            
            //se a thread pertence ao conjunto de multiplos de potencias de dois atual vai realizar o merge
            if(omp_get_thread_num()%((int)pow(2,aux))==0){
                m=r;
                if(omp_get_thread_num()+(int)pow(2,aux)>(omp_get_num_threads()-1))
                    r=size-1;
                else
                    r=( div * (omp_get_thread_num()+(int)pow(2,aux)))-1;
                printf("%d-%d-%d ",omp_get_thread_num(),aux,r);
                merge(V,l,m,r);
                
            }  
            
            #pragma omp barrier
            if (omp_get_thread_num()==0){
                iterate=iterate/2;
                aux+=1;
            }
            #pragma omp barrier
        }
        //merge dos vetores locais
 
    }
    
}
int main(int argc, char **argv)
{
    int* V;//vetor a ser ordenado
    int N;//tmanho do vetor 
    double start,end,cpu_time_used;//contagem do tempo
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