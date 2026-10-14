# Makefile for the CG_analysis dir

CC = g++

CFLAGS = -w -g -std=c++17

OMP_FLAGS = -fopenmp

main: main.cpp analysis.cpp
	$(CC) $(CFLAGS) main.cpp analysis.cpp -o ./bin/main.o

msd: msd.cpp analysis.cpp
	$(CC) $(CFLAGS) msd.cpp analysis.cpp -o ./bin/msd.o

pdf_tavg: pdf_tavg.cpp analysis.cpp
	$(CC) $(CFLAGS) pdf_tavg.cpp analysis.cpp -o ./bin/pdf.o

order: laneOrder.cpp analysis.cpp
	$(CC) $(CFLAGS) laneOrder.cpp analysis.cpp -o ./bin/order.o

rdf1D: rdf1D.cpp analysis.cpp
	$(CC) $(CFLAGS) rdf1D.cpp analysis.cpp -o ./bin/rdf1D.o

rdf2D: rdf2D.cpp analysis.cpp
	$(CC) $(CFLAGS) rdf2D.cpp analysis.cpp -o ./bin/rdf2D.o

disp: displacements.cpp analysis.cpp
	$(CC) $(CFLAGS) displacements.cpp analysis.cpp -o ./bin/disp.o

sk: structure.cpp analysis.cpp
	$(CC) $(CFLAGS) structure.cpp analysis.cpp -o ./bin/sk.o

res: makeRestart.cpp analysis.cpp
	$(CC) $(CFLAGS) makeRestart.cpp analysis.cpp -o ./bin/res.o

merge: mergeTraj.cpp analysis.cpp 
	$(CC) $(CFLAGS) mergeTraj.cpp analysis.cpp -o ./bin/merge.o

potential: computePotentialEnergy.cpp computeInteractions.cpp analysis.cpp ../../LD/LD-cpp/src/interactions.cpp
	$(CC) $(CFLAGS) $(OMP_FLAGS) computePotentialEnergy.cpp computeInteractions.cpp analysis.cpp ../../LD/LD-cpp/src/interactions.cpp -o ./bin/potential.o

temp: temperature.cpp analysis.cpp 
	$(CC) $(CFLAGS) temperature.cpp analysis.cpp -o ./bin/temp.o

power: power.cpp computeInteractions.cpp analysis.cpp ../../LD/LD-cpp/src/interactions.cpp
	$(CC) $(CFLAGS) power.cpp computeInteractions.cpp analysis.cpp ../../LD/LD-cpp/src/interactions.cpp -o ./bin/power.o

collision: collisions.cpp analysis.cpp
	$(CC) $(CFLAGS) collisions.cpp analysis.cpp -o ./bin/coll.o

velocity: velocity.cpp analysis.cpp
	$(CC) $(CFLAGS) velocity.cpp analysis.cpp -o ./bin/velocity.o

vacf: vacf.cpp analysis.cpp
	$(CC) $(CFLAGS) $(OMP_FLAGS) vacf.cpp analysis.cpp -o ./bin/vacf.o

clean:
	rm -f ./bin/*.o
