# options
CFLAGS= -fpic -O3 -pipe -g -c -fopenmp -frounding-math 



all: lib/gpfso.so 



lib/gpfso.so:  functions.o  
	g++ -shared -o lib/gpfso.so  functions.o  -lgsl   -lgslcblas 


functions.o: src/functions.c
	gcc $(CFLAGS) src/functions.c -o functions.o

# clean

clean: 
	rm *.o















