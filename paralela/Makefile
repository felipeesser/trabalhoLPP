all: main.exe
 
main.exe: main.o util.o
	mpicc -o main.exe main.o util.o -lm

main.o: main.c
	mpicc -c main.c -lm

util.o: util.c
	mpicc -c util.c -lm

clean:
	rm -rf *.o main.exe
