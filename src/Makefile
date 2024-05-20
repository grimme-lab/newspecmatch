PROG = ../build/newspecmatch

OBJDIR = ../build
#--------------------------------------------------------------------------
 OSTYPE=LINUXI
 #MKLROOT=/opt/intel/mkl
#--------------------------------------------------------------------------

#modules
MODULES += spectramod.o
MODULES += newspecmatch.o

#mainfile
MAINFILES += main.o 

OBJS = $(MODULES) $(MAINFILES)
vpath % .
OBJS := $(addprefix $(OBJDIR)/, $(OBJS))
#--------------------------------------------------------------------------

ifeq ($(OSTYPE),LINUXI)
  PREOPTS :=
  FC := ifort
  CC := icc
  LINKER = ifort -g -O -static -fopenmp
#  LINKER = ifort -g -O -static -fopenmp  -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
#  LINKER = ifort -static -fopenmp  -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
#  LIBS = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm
#  LIBS := -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a -Wl,--end-group 
  PREFLAG := -E -P
  FFLAGS := -g -O -qopenmp -align array64byte -module $(OBJDIR)
#  FFLAGS := -g -O0
#  FFLAGS := -O -axAVX -qopenmp -align array64byte 
  FFLAGS := -check all -g -debug all -traceback -fpe0
  CCFLAGS := -O -DLINUX
endif

ifeq ($(OSTYPE),LINUXG)
  PREOPTS :=
  FC := gfortran
  CC := gcc
  LINKER := gfortran -static-libgfortran -fopenmp 
  LIBS :=
  PREFLAG :=
  FFLAGS :=  -J$(OBJDIR)
  CCFLAGS := -O -DLINUX
endif


ifeq ($(OSTYPE),SOMETHINGELSE)
  PREOPTS :=
  FC := 
  CC := 
  LINKER := 
  LIBS :=
  PREFLAG := 
  FFLAGS := 
  CCFLAGS := 
endif


.PHONY: all
.PHONY: setup
.PHONY: clean

all: setup
all: $(PROG)
setup: $(OBJDIR)

#--------------------------------------------------------------------------
# Setting the rules to make the object files:
#--------------------------------------------------------------------------
# implizite Regel zur Erzeugung von *.o aus *.F ausschalten
%.o: %.F
fortran.o : $(CUDA)/src/fortran.c
	gcc $(CCFLAGS) -I$(CUDA)/src/ $(INC) -c  $(CUDA)/src/fortran.c

# aus *.F mache ein *.f
%.f: %.F
	@echo "making $@ from $<"
	$(CC) $(PREFLAG) $(PREOPTS) $< -o $@

# aus *.f mache ein *.o
$(OBJDIR)/%.o: %.f
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

$(OBJDIR)/%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

# aus *.c mache ein *.o
$(OBJDIR)/%.o: %.c
	@echo "making $@ from $<"
	$(CC) $(CCFLAGS) -c $< -o $@

# linken
$(PROG): $(OBJS)
	$(LINKER) $(OBJS) $(LIBS) -o $(PROG)

$(OBJDIR):
	mkdir -p $@

#aufraeumen
clean:
	rm -f $(filter %.o,$(OBJS)) $(PROG) 
	rm -f $(OBJDIR)/*.mod *.mod $(PROG) 
	rm -f $(patsubst %.F, %.f, $(wildcard *.F))

