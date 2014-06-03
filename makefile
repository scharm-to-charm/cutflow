# makefile for SUSY ntuple skimmer
# Author: Dan Guest (dguest@cern.ch)
# Created: Wed Jul  4 15:30:40 CEST 2012

# --- set dirs
BIN          := bin
SRC          := src
INC          := include
DICT         := dict
DEP          := $(BIN)


ifndef ROOTCOREDIR
  $(error "couldn't find ROOTCOREDIR")
endif
SUSYTOOLS_INC      := $(shell ./map_libs.sh -i )

#  set search path
vpath %.o    $(BIN)
vpath %.cxx  $(SRC) 
vpath %.hh   $(INC) $(SUSYTOOLS_INC)
vpath %.h    $(INC) $(SUSYTOOLS_INC)
vpath %Dict.h $(DICT)
vpath %Dict.cxx $(DICT)

# --- load in root config
ROOTCFLAGS    := $(shell root-config --cflags)
ROOTLIBS      := $(shell root-config --libs)
ROOTLDFLAGS   := $(shell root-config --ldflags)
ROOTLIBS      += -lTreePlayer   #don't know why this isn't loaded by default
ROOTLIBS      += -lTMVA         #don't know why this isn't loaded by default
ROOTLIBS      += -lXMLParser    #don't know why this isn't loaded by default
ROOTLIBS      += -lEG           #for TParticle


# --- set compiler and flags (roll c options and include paths together)
CXXFLAGS     := -O2 -Wall -fPIC -I$(INC) $(SUSYTOOLS_INC:%=-I%) -g 
COMPILER_NAME := $(notdir ${CXX})
ifeq ($(COMPILER_NAME), g++)
  LDFLAGS      := -Wl,-no-undefined
endif
LIBS         := $(shell ./map_libs.sh -l )

# rootstuff 
CXXFLAGS     += $(ROOTCFLAGS)
LDFLAGS      += $(ROOTLDFLAGS)
LIBS         += $(ROOTLIBS)

# pystuff (roll the linking options and libraries together)
PY_LDFLAGS := $(LDFLAGS)
PY_LDFLAGS += $(PY_LIB)
PY_LDFLAGS += -shared

# dependency flags
DEPINCLUDE := -I$(INC) $(SUSYTOOLS_INC:%=-I%) -I$(shell root-config --incdir)
DEPFLAGS    = -M -MP -MT $(BIN)/$*.o -MT $(DEP)/$*.d $(DEPINCLUDE) $(PY_FLAGS)

# ---- define objects
TOBJ        := SusyBuffer.o 
T_DICTS     := $(TOBJ:.o=Dict.o)
GEN_OBJ     := SmartChain.o CutCounter.o CtagCalibration.o sbottom_functions.o
GEN_OBJ     += PileupReweighting.o
GEN_OBJ     += sbottom_functions.o mctlib.o
EXE_OBJ     := cutflow.o 

ALLDEPOBJ   := $(TOBJ) $(EXE_OBJ) $(GEN_OBJ)
ALLOBJ      := $(ALLDEPOBJ) $(TDICTS)

OUTPUT    := cutflow

all: $(OUTPUT) 

$(OUTPUT): $(ALLOBJ:%=$(BIN)/%)
	@echo "linking $^ --> $@"
	@$(CXX) -o $@ $^ $(LIBS) $(LDFLAGS)

# --------------------------------------------------

# root dictionary generation 
$(DICT)/%Dict.cxx: %.h LinkDef.hh
	@echo making dict $@
	@mkdir -p $(DICT)
	@rm -f $(DICT)/$*Dict.h $(DICT)/$*Dict.cxx 
	@rootcint $@ -c $(INC)/$*.h	
	@sed -i 's,#include "$(INC)/\(.*\)",#include "\1",g' $(DICT)/$*Dict.h

$(BIN)/%Dict.o: $(DICT)/%Dict.cxx 
	@mkdir -p $(BIN)
	@echo compiling dict $@
	@$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -c $< -o $@

# compile rule
$(BIN)/%.o: %.cxx
	@echo compiling $<
	@mkdir -p $(BIN)
	@$(CXX) -c $(CXXFLAGS) $< -o $@

# use auto dependency generation

ifneq ($(MAKECMDGOALS),clean)
  ifneq ($(MAKECMDGOALS),rmdep)
    include $(ALLDEPOBJ:%.o=$(DEP)/%.d)
  endif
endif

$(DEP)/%.d: %.cxx
	@echo making dependencies for $<
	@mkdir -p $(DEP)
	@$(CXX) $(DEPFLAGS) $< -o $@ 

# clean
.PHONY : clean rmdep
CLEANLIST     = *~ *.o *.o~ *.d core 
clean:
	rm -fr $(CLEANLIST) $(CLEANLIST:%=$(BIN)/%) $(CLEANLIST:%=$(DEP)/%)
	rm -fr $(BIN) $(OUTPUT) $(DICT) $(DEP)

rmdep: 
	rm -f $(DEP)/*.d