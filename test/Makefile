HPPATH		= /opt/Healpix/Healpix
HPLIBPATH	= $(HPPATH)/src/cxx/generic_gcc/lib
HPINCPATH	= $(HPPATH)/src/cxx/generic_gcc/include
LAPACKPATH	= /opt/lapack-3.6.0/
MKLPATH		= /opt/intel/compilers_and_libraries/linux/mkl

CC			= g++
CFLAGS		= -fopenmp -O3 -std=c++11
LIBS		= -lcxxsupport -lfftpack -lc_utils -lcurses
HPLIBS		= -lhealpix_cxx -lcfitsio -lsharp
LALIBS		= -llapacke -lcblas -llapack -lrefblas -lgfortran
LALIBPATH	= $(LAPACKPATH)
LAINCPATH	= $(LAPACKPATH)/LAPACKE/include/

MKLINC		= -I$(MKLPATH)/include
MKLLIBS     = -L$(MKLPATH)/lib/intel64 \
			  -L/opt/intel/tbb/lib/intel64/gcc4.7 \
			  -Wl,--start-group \
			  -lmkl_gf_lp64 \
			  -lmkl_core \
			  -lmkl_gnu_thread \
			  -Wl,--end-group 

MINOPT  = -I$(ROOTSYS)/include -L$(ROOTSYS)/lib -lCore -lMinuit 
ROOTOPT = `root-config --ldflags` `root-config --libs` `root-config --cflags`

TARGET2=test_coordtrans
$(TARGET2) : $(TARGET2).cc 
	$(CC) -o $(TARGET2) $(TARGET2).cc  $(CFLAGS) -I$(HPINCPATH) -L$(HPLIBPATH) $(LIBS) $(HPLIBS) $(LIBS)

