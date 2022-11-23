CC = g++ 
PRODUCTFLAGS = -c -std=c++20 -Wall -Wextra -Wpedantic -O3 -march=native -Xpreprocessor -fopenmp

.SUFFIXES: .cpp .o .do

.cpp.o:; $(CC) -o $@ $(PRODUCTFLAGS) $<

coeff_recovery: coefficient_recovery.o rounds_1_to_4.o rounds_5_6.o
	$(CC) -lomp -o coeff_recovery.out $^

coeff_recovery_ubuntu: coefficient_recovery.o rounds_1_to_4.o round_5_6.o
	$(CC) -fopenmp -o superpoly_recovery.out $^

clean:
	rm -f *.o

clean_everything:
	rm -f *.o *.out *.do
