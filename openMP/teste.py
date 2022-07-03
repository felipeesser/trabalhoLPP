import subprocess

for j in range(3):
    for w in range(31):
        total=0.0
        for i in range(3):
            subprocess.call(['./main','../aleatorio/entradas/entrada'+str(w+1),'saida.txt',str(pow(2,j+1))])
            fr = open("tempo_par", "r")
            line=fr.readline().split(' ')
            if (i==0):
                size=line[0]
            total+=float(line[1].rstrip())
            fr.close()
        fw=open("tempo_par_omp_"+str(pow(2,j+1)),"a")
        fw.write(size+' '+str(round(total/3,6))+'\n')
        fw.close()