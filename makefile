FC=gfortran
FCFLAGS = -c -w -fbounds-check

OBJS=	modules.o \
		DIRDFN.o \
        struttura_dati.o \
        DFNcon.o\
        sobol.o\
        halton.o\
		problem.o\
		main_test.o

all: exe

exe: 	$(OBJS)
	$(FC) -o dirdfn $(OBJS)

.SUFFIXES : .f90 .o

.f90.o: $* ; $(FC) $(FCFLAGS) $*.f90

clean:
	rm -f *.o
	rm -f *.mod
	rm -f fort.3
	rm -f dirdfn
