all: main.exe
 
main.exe: main.o util.o
	mpicc -o main.exe main.o util.o

main.o: main.c
	mpicc -c main.c

util.o: util.c
	mpicc -c util.c

clean:
	rm -rf *.o main.exe