# exe name                                                            
#EXE = spaf
EXE = spaf_z

# compiler
CC = gcc
CXX = g++

# preprocessor flags
# 2 domains version
CFLAGS += -DVERSION_Z
CFLAGS += -DTWO_WAY
# incremental pressure correction
# CFLAGS += -DPRE_INCREMENTAL_ROTATIONAL
# use float instead of double in ANN, faster as long as minimum element edge is larger than float machine epsilon
CFLAGS += -DANN_USE_FLOAT

# compiler flags
# CFLAGS += -fopenmp
CFLAGS += -O3 -march=native

# CFLAGS += -g -O0 -Wall -ggdb -g3

# icc optimization AMD64
# CFLAGS += -xW -axW -O3
# CFLAGS += -axP -xW -ipo -O3

# blas library
# CFLAGS += -DUSE_ACCELERATE
# LDFLAGS += -framework Accelerate
# CFLAGS += -DUSE_ACML
# CFLAGS  += -I$(ACML_INCLUDE)
# LDFLAGS += -lacml
# CFLAGS += -DUSE_MKL
# use intel mkl sparse matrix vector product
# CFLAGS += -DUSE_SPARSE_MKL
# LDFLAGS += -lmkl -lguide -lpthread

# taucs library
LDFLAGS += -L$(PWD)/src/taucs -ltaucs

# math library
LDFLAGS += -lm

# source directories
SRC_DIR += src src/nns src/nns/ANN src/misc src/mesh src/solver src/particles src/io src/io/dclog src/fluid src/bc src/linalg src/operators src/sparse src/taucs

# include path
CFLAGS += $(foreach dir,$(SRC_DIR),-I$(dir))  

# misc stuff
CXXFLAGS := $(CFLAGS)
CFLAGS += -std=c99

# OBJ and dependancies directories
OBJ_DIR = .objects
DEP_DIR = .dependencies

# self pointing makefile
CURRENTMAKEFILE=$(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))

#---------------------------------------------------------------------
# Paths where to look
#---------------------------------------------------------------------
vpath %.o $(OBJ_DIR)
vpath %.c $(SRC_DIR)
vpath %.cpp $(SRC_DIR)
vpath %.h $(SRC_DIR)
vpath %.f90 $(SRC_DIR)
vpath %.mod $(OBJ_DIR)
vpath %.d $(DEP_DIR)
#---------------------------------------------------------------------
# Pattern rules
#---------------------------------------------------------------------
%.d: %.cpp $(CURRENTMAKEFILE)
	@echo "Making dependencies for $<"
	@set -e; rm -f $@; \
	$(CXX) -MM $(CXXFLAGS) $< | sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' > $(DEP_DIR)/$@

%.d: %.c $(CURRENTMAKEFILE)
	@echo "Making dependencies for $<"
	@set -e; rm -f $@; \
	$(CC) -MM $(CFLAGS) $< | sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' > $(DEP_DIR)/$@

%.d: %.f90 $(CURRENTMAKEFILE)
	@echo "Making dependencies for $<"
	@set -e; rm -f $@; \
	$(FC) -MM $(FFLAGS) $< | sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' > $(DEP_DIR)/$@

%.o : %.c $(CURRENTMAKEFILE)
	@echo "Compiling $<" ; \
	$(CC) $(CFLAGS) -c $< -o $(OBJ_DIR)/$@

%.o : %.cpp $(CURRENTMAKEFILE)
	@echo "Compiling $<" ; \
	$(CXX) $(CXXFLAGS) -c $< -o $(OBJ_DIR)/$@

%.o : %.f90 $(CURRENTMAKEFILE)
	@echo "Compiling $<" ; \
	$(FC) $(FFLAGS) -c $< -o $(OBJ_DIR)/$@
#---------------------------------------------------------------------
# list of all source files, object files, dependance files, without path
#---------------------------------------------------------------------
SRC = $(foreach dir,$(SRC_DIR),$(patsubst $(dir)/%.c,%.c, $(wildcard $(dir)/*.c) ))
SRC += $(foreach dir,$(SRC_DIR),$(patsubst $(dir)/%.cpp,%.cpp, $(wildcard $(dir)/*.cpp) ))
SRC += $(foreach dir,$(SRC_DIR),$(patsubst $(dir)/%.f90,%.f90, $(wildcard $(dir)/*.f90) ))
FILES = $(basename $(SRC) )
OBJ  = $(addsuffix .o, $(FILES) )
DEPS = $(addsuffix .d, $(FILES) )
#---------------------------------------------------------------------
# Rules                                                               
#---------------------------------------------------------------------
all : directories $(DEPS) $(EXE)

# linking
$(EXE) : $(OBJ)
	@echo "Linking $@" ; \
	cd $(OBJ_DIR) ; \
	$(CXX) $(CXXFLAGS) $(OBJ) $(LDFLAGS) -o ../$(EXE)

# Clean compilation files
.PHONY : clean
clean :
	rm -rf $(OBJ_DIR) $(DEP_DIR)
	rm -f $(EXE) core*

# make sure some directories exist
.PHONY : directories
directories :
	@if [ ! -d $(OBJ_DIR) ] ; then mkdir -p $(OBJ_DIR) ; fi;
	@if [ ! -d $(DEP_DIR) ] ; then mkdir -p $(DEP_DIR) ; fi;

# display flags
.PHONY : cflags
cflags :
	@echo "cflags : $(CFLAGS)"

# display flags
.PHONY : ldflags
ldflags :
	@echo "ldflags : $(LDFLAGS)"

# include the dependencies, except if clean target has been called
ifneq ($(strip $(MAKECMDGOALS)),c)
ifneq ($(strip $(DEPS)),)
-include $(DEPS:%.d=$(DEP_DIR)/%.d)
endif
endif
