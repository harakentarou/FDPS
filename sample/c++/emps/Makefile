PS_PATH = ../../../src/
INC = -I$(PS_PATH)
CC = time g++
CFLAGS = -std=c++17 -O3
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
all:emps.out

emps.out:main.cpp
	$(CC) $(INC) $(CFLAGS) -o $@ $^
