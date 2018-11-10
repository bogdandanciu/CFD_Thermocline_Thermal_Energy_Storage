# Makefile 

CC = g++
CFLAGS = -O3 -Wall 

all: P1 

#P1: P1.o
#	$(CC) -o $@ $@.o $(OBJECTS)  

P1: P1.cpp
	$(CC) $(CFLAGS) P1.cpp -o P1

clean: 
	rm P1 example.txt simData.dat 
