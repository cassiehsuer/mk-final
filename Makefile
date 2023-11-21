CCX = g++
CCXFLAGS = -O3 -funroll-loops -march=native -std=c++11 -pthread -I. -I./include
DEPS = -lntl -lgmp -lfftw3 -lm

all: clean test

clean:
	$(RM) test test.o lwehe.o fft.o sampler.o keygen.o libfinal.a ksk.txt bsk.txt hpk.txt cipher.txt

test: FINAL.h libfinal.a
	$(CCX) $(CCXFLAGS) -o test test.cpp libfinal.a $(DEPS)

libfinal.a: include/params.h lwehe.o keygen.o fft.o sampler.o
	$(AR) -q libfinal.a lwehe.o keygen.o fft.o sampler.o

lwehe.o: include/lwehe.h keygen.o sampler.o src/lwehe.cpp
	$(CCX) $(CCXFLAGS) -c src/lwehe.cpp

keygen.o: include/keygen.h sampler.o fft.o src/keygen.cpp
	$(CCX) $(CCXFLAGS) -c src/keygen.cpp

fft.o: include/fft.h
	$(CCX) $(CCXFLAGS) -c src/fft.cpp

sampler.o: include/sampler.h include/params.h src/sampler.cpp
	$(CCX) $(CCXFLAGS) -c src/sampler.cpp