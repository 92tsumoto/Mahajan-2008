#INCLUDES = syspara.h

bin_PROGRAMS = Mahajan
#Mahajan_CFLAGS = -O2 -Wall -axsse4.2 -mkl -parallel -par-report2 -Winline -inline-level=2 -ip -vec-report2 -opt-mem-layout-trans=3 -pg
Mahajan_CFLAGS = -O2 -Wall -axsse4.2 -mkl -parallel -par-report2 -vec-report2 -pg

Mahajan_INCLUDE = -I./include
Mahajan_LDADD = -L/opt/intel/mkl/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core ./lib/xhplot.o

Mahajan_SOURCES = ./src/mahajan08.c ./src/var_consts.c ./src/mem.c\
	./src/win_out.c ./src/eular.c ./src/input_para.c ./src/system.c\
	./src/eventloop.c ./src/exptable.c ./src/data_out.c ./src/syspara.h

