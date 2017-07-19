
MAKEFILE_OPT = $(PWD)/Makefile_opt.in
include $(MAKEFILE_OPT)

MAKEFILE_H5  = $(PWD)/Makefile_dir.in
include $(MAKEFILE_H5)

TEMPLATES = bexp bx3d earth fieldloop flock flock_grmhd isentropic jupiter kepler kh mri2 rotor shear shocktube spinring spread vortex sorathia_grmhd blast_grmhd fieldloop_grmhd bl kep_ring

GIT_VERSION = $(shell git describe --dirty --always --tags)

OPT_DEFS = -DGIT_VERSION=\"$(GIT_VERSION)\"
OPT_DEFS += -DINITIAL=\"$(INITIAL)\"
OPT_DEFS += -DHYDRO=\"$(HYDRO)\"
OPT_DEFS += -DBOUNDARY=\"$(BOUNDARY)\"
OPT_DEFS += -DOUTPUT=\"$(OUTPUT)\"
OPT_DEFS += -DRESTART=\"$(RESTART)\"
OPT_DEFS += -DPLANETS=\"$(PLANETS)\"
OPT_DEFS += -DHLLD=\"$(HLLD)\"
OPT_DEFS += -DANALYSIS=\"$(ANALYSIS)\"
OPT_DEFS += -DMETRIC=\"$(METRIC)\"
OPT_DEFS += -DFRAME=\"$(FRAME)\"
OPT_DEFS += -DNUM_C=$(NUM_C)
OPT_DEFS += -DNUM_N=$(NUM_N)
OPT_DEFS += -DCT_MODE=$(CT_MODE)

FLAGS = -O3 -Wall -g $(OPT_DEFS)

INC = -I$(H55)/include
LIB = -L$(H55)/lib -lhdf5 -lm

OBJ = main.o readpar.o timestep.o onestep.o riemann.o mpisetup.o gridsetup.o domain.o misc.o geometry.o faces_alt.o exchange.o plm.o report.o profiler.o planet.o omega.o analysis.o bfields.o $(HLLD).o rotframe.o boundary_functions.o $(INITIAL).o $(OUTPUT).o $(HYDRO).o $(BOUNDARY).o $(RESTART).o $(PLANETS).o $(METRIC).o $(FRAME).o calc.o $(ANALYSIS).o  noise.o #snapshot.o

default: disco

.PHONY: $(TEMPLATES)

$(TEMPLATES):
	cp Templates/$@.par in.par
	cp Templates/$@.in Makefile_opt.in
	make clean
	make

%.o: %.c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c $<

calc.o: Calc/*.c Calc/calc.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Calc/*.c -o calc.o

$(TIMESTEP).o: Timestep/$(TIMESTEP).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Timestep/$(TIMESTEP).c

$(INITIAL).o : Initial/$(INITIAL).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Initial/$(INITIAL).c

$(HYDRO).o : Hydro/$(HYDRO).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Hydro/$(HYDRO).c

$(PLANETS).o : Planets/$(PLANETS).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Planets/$(PLANETS).c

$(BOUNDARY).o : Boundary/$(BOUNDARY).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Boundary/$(BOUNDARY).c

$(OUTPUT).o : Output/$(OUTPUT).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Output/$(OUTPUT).c

$(RESTART).o : Restart/$(RESTART).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Restart/$(RESTART).c

$(ANALYSIS).o : Diagnostics/$(ANALYSIS).c paul.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Diagnostics/$(ANALYSIS).c

$(METRIC).o : Hydro/Metric/$(METRIC).c paul.h Hydro/metric.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Hydro/Metric/$(METRIC).c

$(FRAME).o : Hydro/Frame/$(FRAME).c paul.h Hydro/frame.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c Hydro/Frame/$(FRAME).c

disco: $(OBJ) paul.h
	$(CC) $(FLAGS) -o disco $(OBJ) $(LOCAL_LDFLAGS) $(LIB)

clean:
	rm -f *.o disco
