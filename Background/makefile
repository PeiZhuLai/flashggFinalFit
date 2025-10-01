# name of the library
LIBNAME = BackgroundProfileFitting

# Necessary to use shell built-in commands
SHELL = bash

# figure out os
UNAME   := $(shell uname)
BASEDIR := $(shell pwd)

# create directories
$(shell mkdir -p bin)
$(shell mkdir -p lib)
$(shell mkdir -p obj)

# ---- turn off built-in implicit rules to avoid accidental link steps ----
.SUFFIXES:
MAKEFLAGS += --no-builtin-rules

# includes & libs
USERINCLUDES += -I$(ROOTSYS)/include
# also bring in root-config cflags (includes -I and necessary macros)
USERINCLUDES += $(shell root-config --cflags)

ifeq ($(UNAME),Darwin)
  $(info Found OS Darwin - using mac osx build)
  BOOST_PATH=/Applications/boost_1_49_0
  USERINCLUDES += -I$(BOOST_PATH)/include/
  USERLIBS += -L$(BOOST_PATH)/libs/ -lboost_regex -lboost_program_options
endif

ifeq ($(UNAME),Linux)
  $(info Found OS Linux - using linux build)
  BOOST_PATH=$(shell scram tool info boost | awk -F= '/BOOST_BASE/ {  print $$2  }')
  USERINCLUDES += -I$(BOOST_PATH)/include/
  USERINCLUDES += -I $(CMSSW_BASE)/src/ -I $(CMSSW_RELEASE_BASE)/src/
  USERLIBS += -L$(BOOST_PATH)/lib/ -lboost_regex -lboost_program_options
  USERINCLUDES += -I$(ROOFITSYS)/include
  USERLIBS += -L$(ROOFITSYS)/lib -lRooFit -lRooFitCore -lRooStats
  USERLIBS += -L$(CMSSW_BASE)/lib
endif

USERLIBS += $(shell root-config --glibs) -lTMVA -lMinuit -lFoam -lRooFit -lRooFitCore -lMathMore -lRooStats
CXXFLAGS = -Wall -g -Wl,-rpath=$(BASEDIR)/lib -std=c++17
# global hardening: all compiles with -fPIC and disable LTO to avoid mixed-link issues
CXXFLAGS += -fPIC -fno-lto
LDFLAGS = -shared -W

CXX = g++
LD  = g++

# dirs
LIBDIR  = $(BASEDIR)/lib
EXEDIR  = $(BASEDIR)/bin
MACRODIR= $(BASEDIR)/src
SRCDIR  = $(BASEDIR)/src
OBJDIR  = $(BASEDIR)/obj
INCDIR  = $(BASEDIR)/interface
TESTDIR = $(BASEDIR)/test
DOCDIR  = $(BASEDIR)/docs

OBJ_EXT  = o
TEST_EXT = cpp

# include paths
INCLUDE  += -I $(INCDIR)
INCLUDE  += $(USERINCLUDES)
CXXFLAGS += $(INCLUDE)

LIBS += $(USERLIBS)
LIBS += -L$(CMSSW_BASE)/lib/$(SCRAM_ARCH)/ -lHiggsAnalysisCombinedLimit

# files
ROOT_DICT = $(OBJDIR)/RootDict.cxx
ROOT_OBJ  = $(OBJDIR)/RootDict.$(OBJ_EXT)

SRCS    := $(wildcard $(SRCDIR)/*.cc)
EXES    := $(wildcard $(TESTDIR)/*.$(TEST_EXT))
HEADERS := $(wildcard $(INCDIR)/*.h)

# robust mappings (avoid whitespace issues of subst)
OBJS := $(patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.$(OBJ_EXT),$(SRCS))
BINS := $(patsubst $(TESTDIR)/%.$(TEST_EXT),$(EXEDIR)/%,$(EXES))

# guards: forbid treating these as link targets
src/%: ;
interface/%: ;
obj/%: ;

.PHONY: all clean lib vars

all: $(BINS)

$(EXEDIR)/%: $(TESTDIR)/%.$(TEST_EXT) $(LIBDIR)/lib$(LIBNAME).so $(INCDIR)/*.h
	@echo Compiling executable $@
	@$(CXX) -o $@ $(CXXFLAGS) $< -L$(LIBDIR) -l$(LIBNAME) $(LIBS)

$(OBJDIR)/%.$(OBJ_EXT): $(SRCDIR)/%.cc $(ROOT_DICT)
	@echo Making object $@
	@$(CXX) $(CXXFLAGS) -c $< -o $@

$(LIBDIR)/lib$(LIBNAME).so: $(OBJS) $(ROOT_OBJ)
	@echo Building shared library $@
	@$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

# 1) 由頭文件產生 RootDict.cxx
$(ROOT_DICT): $(HEADERS)
	@echo Making dictionary $@
	@rootcint -f $@ -c -I$(INCDIR) $(HEADERS)

# 2) 由 RootDict.cxx 編譯出 RootDict.o
$(ROOT_OBJ): $(ROOT_DICT)
	@echo Compiling dictionary object $@
	@$(CXX) $(CXXFLAGS) -fPIC -c $(ROOT_DICT) -o $@
	
lib: $(LIBDIR)/lib$(LIBNAME).so

vars:
	@echo "LIBS:            " $(LIBS)
	@echo "CXXFLAGS:        " $(CXXFLAGS)
	@echo "Header files:    " $(HEADERS)
	@echo "Source files:    " $(SRCS)
	@echo "Object files:    " $(OBJS)
	@echo "Executables:     " $(BINS)
	@echo "ROOT_DICT:       " $(ROOT_DICT)
	@echo "ROOT_OBJ:        " $(ROOT_OBJ)

clean:
	@echo "Cleaning:"
	@echo "  " $(OBJS)
	@echo "  " $(LIBDIR)/lib$(LIBNAME).so
	@echo "  " $(BINS)
	@echo "  " $(ROOT_DICT)
	@rm -rf $(OBJS) $(LIBDIR)/lib$(LIBNAME).so $(BINS) $(ROOT_DICT) $(ROOT_OBJ) $(OBJDIR)/RootDict_rdict.pcm
