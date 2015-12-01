# ======================================
# Template C++11 Makefile
# author: Gustavo Peixoto de Oliveira
# e-mail: gustavo.oliveira@uerj.br
# ======================================

# compiler
CXX = clang++
#CXX = g++

# flags
CXXFLAGS = -g -std=c++11 -Wall 
LFLAGS = -lm

# directories
BASE_DIR = .

# inclusions
INCLUDES += -I${BASE_DIR} 
INCLUDES += -I${MATLAB_DIR}/extern/include
INCLUDES += -I${SNAP_DIR}/snap-core
INCLUDES += -I${SNAP_DIR}/snap-exp
INCLUDES += -I${SNAP_DIR}/snap-adv
INCLUDES += -I${SNAP_DIR}/glib-core
INCLUDES += -I${SNAP_DIR}/glib-adv

UNAME = $(shell uname -s)

# Linux
ifeq ($(findstring Linux,${UNAME})),Linux)
 CINCLUDE = -I$(MATLAB_DIR)/extern/include -I. -shared
 EXT = mexa64
endif

# Mac OS
ifeq ($(findstring Darwin,${UNAME})),Darwin)
 CINCLUDE = -L$(MATLAB_DIR)/bin/maci64 -I. -shared -lmx -lmex -lmat
 EXT = mexmaci64
 CXXFLAGS += -std=c++11
endif

# get list of all .cpp files in the current dir
SRC = $(wildcard ${BASE_DIR}/*.cpp) 
SRCH = $(wildcard ${BASE_DIR}/*.h) 

# substitution reference (create .o objects): 
OBJ = $(SRC:%.cpp=%.o)
OBJ2 = $(SNAP_DIR)/snap-core/Snap.o

#=== rules 

all: graphMetrics
	@echo \
	 Executable compiled.

# === right
%.o: %.cpp ${SRCH}
	${CXX} ${CXXFLAGS} ${INCLUDES} -c $< -o $@ 

# executables
graphMetrics: ./main/mainGraphMetrics.o ${OBJ} ${OBJ2}
	${CXX} ${OBJ} ${OBJ2} ${CXXFLAGS} ${LFLAGS} $< -o $@

#=== phony variables

.PHONY: clean
#clean: 

# clean everything 
.PHONY: deepclean
deepclean:
	@echo
	@rm -f exec
	@rm -f ${OBJ} 
	@find ${BASE_DIR} -name "?*.o" -exec rm {} \; 


# clean post-processing dir
.PHONY: cleanPP
cleanPP:
	@find ${PP_DIR} -type f -name "*.*" -exec rm -f {} \;
	@echo Post-processing dir cleaned.

# clean documentation dir
.PHONY: cleandoc
cleandoc:
	@rm -fr ${DOCS_DIR} 
	@echo Documentation dir cleaned.


#=== help 
# automatic variables
# $@ - the file name of the target of the rule
# $< - the name of the first prerequisite (here, source files)
#

# Makefile debugging tools
# - warnings:
# $(warning This is a message that can be placed everywhere in makefile) 
# - prints:
#   'make --just-print' or 'make -n'  
#   'make -p' (data-base)  
#   'make --warn-undefined-variables'
#   'make --debug'
