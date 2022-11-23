CC = g++
PRODUCTFLAGS = -c -std=c++20 -Wall -Wextra -Wpedantic -O3 -march=native -Xpreprocessor -fopenmp 


.SUFFIXES: .cpp .o .do

.cpp.o:; $(CC) -o $@ $(PRODUCTFLAGS) $<

values_recovery: values_recovery.o cube_sum.o permutation.o
	$(CC) -lomp -o values_recovery.out $^

values_recovery_ubuntu: values_recovery.o cube_sum.o permutation.o
	$(CC) -fopenmp -o values_recovery.out $^

# Clean deletes .o files, clean_everything cleans everything, obviously
clean:
	rm -f  *.o
clean_everything:
	rm -f *.o *.out *.do
