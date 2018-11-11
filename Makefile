CC=g++
OBJ = main.o functions.o
HEADER = include/functions.h
CFLAGS = -c -Wall -O3 
CPPFLAGS = -Iinclude 
.PHONY: all clean

run_main: $(OBJ)
	$(CC) $(OBJ) -o $@

main.o: src/main.cpp $(HEADER)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -o $@

functions.o: src/functions.cpp $(HEADER)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -o $@

clean:
	rm *.o run_main

clean_all:
	rm *.dat




























