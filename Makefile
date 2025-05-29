# Makefile for compiling asohf.x with gfortran or ifort
# Usage: make COMP=1 [READER=2]

# Source files
SRC = particles.f asohf.f
OUT = asohf.x

$(info COMP=$(COMP))

# Default READER value: MASCLET
READER := $(if $(READER),$(READER),1)
$(info READER=$(READER))

HDF5 = 0 
ifeq ($(READER),3)
  HDF5=1
endif
$(info HDF5=$(HDF5))

# Compiler selection
ifeq ($(COMP),1)
    FC = gfortran
    FLAGS = -O3 -march=native -fopenmp -mcmodel=medium -funroll-all-loops -fprefetch-loop-arrays -mieee-fp -ftree-vectorize -cpp
else ifeq ($(COMP),2)
    FC = ifort
    FLAGS = -O3 -mcmodel=medium -qopenmp -shared-intel -fp-model consistent -ipo -xHost -cpp
endif

# Preprocessor flags
DEFINES = -Dreader=$(READER)

# Add HDF5 if READER=2
ifeq ($(HDF5),1)
	INC += -I/usr/include/hdf5/serial/
	LIBS += -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5_serial_fortran -lhdf5hl_fortran
endif

# Default target
all:
	$(FC) $(FLAGS) $(DEFINES) $(SRC) $(INC) -o $(OUT) $(LIBS)

# Clean
clean:
	rm -f $(OUT) particles.mod

info:
	@echo ""
	@echo "***********************************************************"
	@echo "***                        ASOHF                        ***"
	@echo "***********************************************************"
	@echo "* David Vallés-Pérez et al., Universitat de València, 2022*"
	@echo "***********************************************************"
	@echo "*** FLAGS for compiling ***"
	@echo "- COMP (mandatory): chooses the compiler options; provided "
	@echo "                    are two examples, but you might need " 
	@echo "                    to modify them for your own machine(s)."
	@echo "   In particular, you might need to change the paths to the "
	@echo "     hdf5 library, if you are using it (e.g. READER=3)."
	@echo ""
	@echo "- READER (optional, default: 1): which code to read the"
	@echo "                    outputs from. Options:"
	@echo "       0: General reader (see format specification in docs)"
	@echo "       1: MASCLET reader (default)"
	@echo "       2: Gadget unformatted reader"
	@echo "       3: Gizmo/Gadget/Arepo HDF5 reader"
	@echo ""
