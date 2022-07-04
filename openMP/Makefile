all: main.exe
 
main.exe: main.o util.o
	gcc -o main.exe main.o util.o -O3 -lm -fopenmp

main.o: main.c
	gcc -c main.c -O3 -lm -fopenmp

util.o: util.c
	gcc -c util.c -O3 -lm -fopenmp

clean:
	rm -rf *.o main.exe