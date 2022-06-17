#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
void geraArq(char *nome_arq,int tam){
    FILE* arq = fopen(nome_arq, "w");
    fprintf(arq,"%d\n",tam);
    srand(time(0));
    for (int i = 0; i < tam; ++i)
    {
       if (i<tam-1)
           fprintf(arq,"%d\n",rand()%(10*tam));
       else
            fprintf(arq,"%d",rand()%(10*tam));
       
    }
    fclose(arq);
}
int main(int argc, char **argv){
   char nome[20]="./entradas/entrada";
   char num_trial[50];
   char aux[70];
   int tam_inicio= atoi(argv[1]);
   int tam_passo=atoi(argv[2]);
   int tam=tam_inicio;
   for (int i = 1; i <= 31; i++)
   {
       strcpy(aux,nome);
       sprintf(num_trial,"%d",i);
       geraArq(strcat(aux,num_trial),tam);
       tam+=tam_passo;
   }
   
   
}