#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "util.h"

void merge(int arr[], int l, int m, int r)
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
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
  
    /* Copy the remaining elements of L[], if there
    are any */
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }
  
    /* Copy the remaining elements of R[], if there
    are any */
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
    free(L);
    free(R);
}
  
void mergeSort(int arr[], int l, int r)
{
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l + (r - l) / 2;
  
        // Sort first and second halves
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);
  
        merge(arr, l, m, r);
    }
}
  
void printArray(int A[], int size)
{
    int i;
    for (i = 0; i < size; i++)
        printf("%d ", A[i]);
    printf("\n");
}

int main(int argc, char **argv)
{
    FILE* f;
    int* V;//vetor a ser ordenado
    int N;//tmanho do vetor 
    double start,end,cpu_time_used;//contagem do tempo
    char *arq_entrada = argv[1], 
    *arq_saida = argv[2];
    ler_entrada(arq_entrada,&V,&N);
    start=clock();
    mergeSort(V, 0,N - 1);
    end=clock();
    cpu_time_used=((double) (end - start)) / CLOCKS_PER_SEC;
    escrever_saida(arq_saida,V,N,cpu_time_used);
    liberar(&V);
    return 0;
}