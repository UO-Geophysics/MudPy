#!/bin/bash
gfortran -O  -c -o fk.o fk.f
gfortran -O  -c -o kernel.o kernel.f
gfortran -O  -c -o prop.o prop.f
gfortran -O  -c -o source.o source.f
cpp -traditional-cpp bessel.FF > bessel.f
gfortran -O  -c -o bessel.o bessel.f
gcc -O   -c -o fft.o fft.c
gcc -O   -c -o Complex.o Complex.c
gcc -O   -c -o sacio.o sacio.c
gfortran -O  -c -o haskell.o haskell.f
gfortran -O   -o fk fk.o kernel.o prop.o source.o bessel.o fft.o Complex.o sacio.o haskell.o -lm
gcc -O   -c -o syn.o syn.c
gcc -O   -c -o radiats.o radiats.c
gfortran -O  -c -o futterman.o futterman.f
gfortran -O   -o syn syn.o fft.o Complex.o sacio.o radiats.o futterman.o  -lm
gfortran -O  -c -o st_haskell.o st_haskell.f
gfortran -O   -o st_fk fk.o kernel.o prop.o source.o bessel.o fft.o Complex.o sacio.o st_haskell.o -lm
gfortran -O  -c -o trav.o trav.f
gfortran -O  -c -o tau_p.o tau_p.f
gfortran -O   -o trav trav.o tau_p.o -lm
gcc -O   -c -o sachd.o sachd.c
gcc -O    -o sachd sachd.o sacio.o -lm
