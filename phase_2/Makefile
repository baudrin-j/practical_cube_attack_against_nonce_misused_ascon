CC = g++ 
PRODUCTFLAGS = -c -std=c++20 -Wall -Wextra -Wpedantic -O3 -march=native -Xpreprocessor -fopenmp

.SUFFIXES: .cpp .o .do

.cpp.o:; $(CC) -o $@ $(PRODUCTFLAGS) $<

phase_2: coefficient_recovery/coefficient_recovery.o coefficient_recovery/rounds_1_to_4.o coefficient_recovery/rounds_5_6.o values_recovery/permutation.o values_recovery/cube_sum.o values_recovery/values_recovery.o
	$(CC) -lomp -o phase_2.out $^

phase_2_ubuntu:coefficient_recovery/coefficient_recovery.o coefficient_recovery/rounds_1_to_4.o coefficient_recovery/rounds_5_6.o values_recovery/permutation.o values_recovery/cube_sum.o values_recovery/values_recovery.o
	$(CC) -fopenmp -o phase_2.out $^

clean:
	find . -name '*.o' -delete

clean_everything:
	find . -name '*.swp' -delete
	rm -f *.out
