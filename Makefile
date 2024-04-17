#
# Makefile basis
#

ifndef PROG
# PROG     = $(shell grep -i program *.f | awk '$$2=="program"&&$$4==""{print $$1}' | sed s/.f://)
 NFILES   = $(shell ls -1 *.f | wc -l)
ifeq ($(NFILES),1)
 PROG     = $(shell grep -i program *.f | awk '$$1=="program"&&$$3==""{print $$2}' )
else
 PROG     = $(shell grep -i program *.f | awk '$$2=="program"&&$$4==""{print $$3}' )
endif
endif
 SOURCES  = $(wildcard *.[f])
 INC      = $(wildcard *.[inc])
 OBJ      = $(shell echo "$(SOURCES)" | sed s/f\ /o\ /g)

 ARCH = $(shell uname -m)

ifdef DEBUG
# EXE   = $(PROG)_$(OSNAME)_$(ARCH)_debug
 EXE   = $(PROG)_$(ARCH)_debug
else
# EXE   = $(PROG)_$(OSNAME)_$(ARCH)
 EXE   = $(PROG)_$(ARCH)
endif

# Linux compiler
 COMPILER = gfortran

# Libraries
ifndef SYSLIBS
 SYSLIBS  = -lnsl
endif

ifndef CERNLIBS
 CERNLIBS = -lpdflib804 -lmathlib -lphtools -lgeant321 \
 -lpacklib -lkernlib -lpawlib
endif

ifndef CLASLIBS
 CLASLIBS = -lbankdefs -lc_bos_io -lfputil \
 -lrecutl -lbos -lfpack -lcc -lmapmanager -lclasutil
endif

# LIBS = $(SYSLIBS) -L$(CLAS_LIB) $(CLASLIBS) -L$(CERNLIB) $(CERNLIBS)
# use only static CERNLIBS (otherwise gives error) and CLASLIBS
# LIBS = $(SYSLIBS) -Wl,-Bstatic -L$(CLAS_LIB) $(CLASLIBS) -L$(CERNLIB) $(CERNLIBS) -Wl,-Bdynamic
# at Ubuntu 22 we have not CERNLIB or CLASLIB so far
LIBS = $(SYSLIBS)

ifndef F77OPT
ifeq "$(ARCH)" "x86_64"
ifdef DEBUG
 F77OPT = -O2 -march=native -mtune=native -m64 -DLinux -fno-automatic -fbounds-check \
 -funroll-all-loops \
 -fdollar-ok -ffixed-line-length-none -fno-second-underscore \
 -Wunused -Wuninitialized -g
else
# F77OPT = -O2 -march=core2 -mtune=core2 -m64 -DLinux -fno-automatic -fbounds-check
# F77OPT = -O2 -march=native -mtune=native -m64 -DLinux -fno-automatic -fbounds-check 
 F77OPT = -O2 -march=native -mtune=native -m64 -DLinux -fno-automatic -fbounds-check \
 -funroll-all-loops \
 -fdollar-ok -ffixed-line-length-none -fno-second-underscore \
 -Wunused -Wuninitialized
endif
else
ifdef DEBUG
 F77OPT = -O2 -march=i686 -m32 -DLinux -fno-automatic -fbounds-check -fno-f2c \
 -funroll-all-loops \
 -fdollar-ok -ffixed-line-length-none -fno-second-underscore \
 -Wunused -Wuninitialized -g
else
 F77OPT = -O2 -march=i686 -m32 -DLinux -fno-automatic -fbounds-check -fno-f2c \
 -funroll-all-loops \
 -fdollar-ok -ffixed-line-length-none -fno-second-underscore \
 -Wunused -Wuninitialized
endif
endif
endif

%.o:    %.f $(INC)
	$(COMPILER) $(F77OPT) $(CPPOPT) -c $< -o $@

all:    $(EXE)

$(EXE): $(OBJ) $(INC)
	$(COMPILER) $(F77OPT) -o $(EXE) $(OBJ) $(LIBS)

clean:
	rm -f  *.o core

delete:
	rm -f $(EXE)

