import subprocess

for w in range(31):
    total=0.0
    for i in range(3):
        subprocess.call(['./merge','../aleatorio/entradas/entrada'+str(w+1),'saida.txt'])
        fr = open("tempo_seq", "r")
        line=fr.readline().split(' ')
        if (i==0):
            size=line[0]
        total+=float(line[1].rstrip())
        fr.close()
    fw=open("tempo_seq_1","a")
    fw.write(size+' '+str(round(total/3,6))+'\n')
    fw.close()