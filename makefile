comp:=gfortran
pattern:=*.f95
paral:=-fopenmp
opt:=-O2
source:=$(wildcard $(pattern))
   obj:=$(patsubst %.f95, %.o, $(source))
prog: $(obj)
	$(comp) $(opt) $^ -o $@
	rm -f *.o
%.o %.mod : %.f95
	$(comp) -c $(opt) $<
main.o : work_program.mod work_function.mod
work_program.o : work_function.mod
clear:
	rm -f *.o *.mod prog
