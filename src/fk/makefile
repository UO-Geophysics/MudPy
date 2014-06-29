optimize=-O
FFLAGS=$(optimize) -m64
CFLAGS=$(optimize) -m64
#if SAC library has been installed, uncomment the next two lines
#CFLAGS=$(optimize) -DSAC_LIB -m32
#SACLIB=-L/Users/dmelgarm/bin/sac/lib -lsac -lsacio
#

CC=gcc
FC=gfortran

SUBS = fft.o Complex.o sacio.o
FKSUBS = fk.o kernel.o prop.o source.o bessel.o $(SUBS)

all: fk syn st_fk trav

syn: syn.o ${SUBS} radiats.o futterman.o
	${LINK.f} -o $@ $^ ${SACLIB} -lm

fk: ${FKSUBS} haskell.o
	${LINK.f} -o $@ $^

st_fk: ${FKSUBS} st_haskell.o
	${LINK.f} -o $@ $^

trav: trav.o tau_p.o
	$(LINK.f) -o $@ trav.o tau_p.o

bessel.f: bessel.FF
	cpp -traditional-cpp $< > $@

clean:
	rm -f *.o bessel.f
