# ======================================================================
# NOTE: make all
# -----
#
# ======================================================================
.PHONY: prep clean vars

ROOTCINT      = $(ROOTSYS)/bin/rootcint
ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

ifdef SCRAM_ARCH
  CXX         := $(shell scramv1 tool info gcc-cxxcompiler | grep CXX= | sed s/CXX=//)
else
  ifdef CXXCOMPILER
    CXX         := $(CXXCOMPILER)
  else
    CXX         := c++
  endif
endif

ROOT6 = $(shell $(ROOTSYS)/bin/root-config --version)
ifneq ($(findstring 6.,$(ROOT6)),)
ROOT6         = "yes"
else
ROOT6         = "no"
endif

CXXFLAGS      = -g -O0 -Wall -fPIC -pipe
LD            = $(CXX)
LDFLAGS       = -g
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
GLIBS         = $(filter-out -lz, $(ROOTGLIBS)) -lTMVA -lTMVAGui -lRooFitCore -lRooFit
EXTHEADERS    = -I..


# -- Files
#CSTUFF =
READER = treeReader01.o

DICT = ${READER:.o=Dict.o}
DICTHEADERS = ${READER:.o=Dict.h}
DICTSOURCES = ${READER:.o=Dict.cc}

# -- Default rules
$(addprefix obj/,%.o) : %.cc %.hh %.icc
	$(CXX) $(CXXFLAGS) $(EXTHEADERS) -c $< -o $@

$(addprefix obj/,%.o) : %.cc %.hh
	$(CXX) $(CXXFLAGS) $(EXTHEADERS) -c $< -o $@

$(addprefix obj/,%.o) : %.cc
	$(CXX) $(CXXFLAGS) $(EXTHEADERS) -c $< -o $@

%Dict.cc : %.hh %LinkDef.h
	$(ROOTCINT) -f $@ -c $(EXTHEADERS) $^
ifeq ($(ROOT6),"yes")
	mv $(subst Dict.cc,Dict_rdict.pcm,$@) lib
endif

%Dict.cc : %.hh
	$(ROOTCINT) -f $@ -c $(EXTHEADERS) $<
ifeq ($(ROOT6),"yes")
	mv $(subst Dict.cc,Dict_rdict.pcm,$@) lib
endif

$(addprefix lib/,%.pcm) :
	cd lib && ln -s $< && cd -


# -- Targets
all: vars prep lib bin

lib: $(addprefix obj/,$(CSTUFF) $(READER)  $(DICT))
	$(CXX) $(SOFLAGS) $(addprefix obj/,$(CSTUFF) $(READER) $(DICT)) -o lib/libAnaClasses.so $(GLIBS) -lMinuit ../build/libClassesDict.so

bin: lib obj/runTreeReader01.o
	$(LD) $(LDFLAGS) -o bin/runTreeReader01 $(GLIBS) obj/runTreeReader01.o lib/libAnaClasses.so ../build/libClassesDict.so

# -- preparatory setup
prep:
	mkdir -p obj bin lib
	cd lib && ln -f -s ../build/libClassesDict.so && cd -

# -- clean up
clean:
	rm -f $(addprefix obj/,$(READER) $(DICT) $(CSTUFF) runTreeReader01.o)
	rm -f $(DICTHEADERS)  $(DICTSOURCES)
	rm -f bin/*
	rm -f lib/*

cleanresults:
	rm -f results/*.root results/*.pdf results/*.tex results/*.tex.old
