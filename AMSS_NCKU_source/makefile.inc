
## filein  = -I/usr/include -I/usr/lib/x86_64-linux-gnu/mpich/include -I/usr/lib/x86_64-linux-gnu/openmpi/lib/

# Use the compiler’s default C++ include paths to avoid mixing libstdc++ versions.
# Hard-coding /usr/include/c++/<version> caused GCC 13 to see headers from two
# different releases (type_traits from 11 with c++config from 13), which broke
# the _GLIBCXX20_DEPRECATED macro expansion. Only add the extra project/library
# headers that are actually required.

# Let the linker pick the matching libgfortran from the active toolchain to
# avoid pinning to a specific GCC major version.
LDLIBS  =  -lgfortran 

#CXXAPPFLAGS =  -Og -Wno-deprecated -Dfortran3 -Dnewc -fopenmp 
#f90appflags = -Og -cpp -ffast-math -march=native -mtune=native
CXXAPPFLAGS = -Ofast -Wno-deprecated -Dfortran3 -Dnewc -fopenmp -ffast-math -funroll-loops -flto -fno-math-errno -ffp-contract=fast -fno-fat-lto-objects -fipa-pta
f90appflags = -O3 -x f95-cpp-input -fopenmp -ffree-line-length-none -ffast-math -march=native -mtune=native -funroll-loops -flto -fno-math-errno -ffp-contract=fast -fno-fat-lto-objects 
f90          = gfortran 
f77          = gfortran 
CXX          = mpic++
CC           = gcc
CLINKER      = mpic++ 

Cu = nvcc
CUDA_LIB_PATH = -L/usr/lib/cuda/lib64 -I/usr/include -I/usr/lib/cuda/include
#CUDA_APP_FLAGS = -c -g -O3 --ptxas-options=-v -arch compute_13 -code compute_13,sm_13 -Dfortran3 -Dnewc
CUDA_APP_FLAGS = -c -g -O3 --ptxas-options=-v -Dfortran3 -Dnewc