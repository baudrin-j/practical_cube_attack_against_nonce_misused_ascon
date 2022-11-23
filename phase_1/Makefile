CC = g++
PRODUCTFLAGS = -c -std=c++17 -Wall -Wextra -Wpedantic -O3 -march=native -Xpreprocessor -fopenmp 

.SUFFIXES: .cpp .o

.cpp.o:; $(CC) -o $@ $(PRODUCTFLAGS) $<

phase_1_verif: phase_1_verification.o random.o cube_sum.o permutation.o
	$(CC) -lomp -o phase_1_verif.out $^

phase_1_verif_ubuntu: phase_1_verification.o random.o cube_sum.o permutation.o
	$(CC) -fopenmp -o phase_1_verif.out $^

# Clean deletes .o files, clean_everything cleans everything, obviously
clean:
	rm -f  *.o
clean_everything:
	rm -f *.o *.out
