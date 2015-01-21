# name of the library
LIBNAME = CMSDASAnalysis

#Necessary to use shell built-in commands
SHELL=bash

# Define include paths
USERINCLUDES += -I$(ROOTSYS)/include/
USERINCLUDES += -Iinterface/
USERINCLUDES += -isystem $(BOOSTSYS)/include/

# Define libraries to link
USERLIBS += $(shell root-config --glibs) -lTMVA
#USERLIBS += -L$(BOOSTSYS)/lib/ -lboost_regex -lboost_program_options -lboost_filesyestem 

#CXXFLAGS = -Wall -W -Wno-unused-function -Wno-parentheses -Wno-char-subscripts -Wno-unused-parameter -O2 
CXXFLAGS = -Wl,--no-as-needed $(shell root-config --glibs) -Wall -W -O2 -std=c++0x  # -std=c++11 
#CXXFLAGS = -Wall -W -O2 -std=c++11 
LDFLAGS = -shared -Wall -W


# If possible we'll use the clang compiler, it's faster and gives more helpful error messages
# If it's not available, then fallback to gcc.
CXX=g++
LD=g++
# CLANGPATH := $(shell type -p clang++)
# ifeq ($(CLANGPATH),)
# $(warning clang++ not found, reverting to g++!)
# CXX=g++
# LD=g++
# endif

CXXFLAGS += $(USERINCLUDES)
LIBS += $(USERLIBS)

# A list of directories
BASEDIR = $(shell pwd)
LIBDIR = $(BASEDIR)/lib
EXEDIR = $(BASEDIR)/bin
SRCDIR = $(BASEDIR)/src
OBJDIR = $(BASEDIR)/obj
TESTDIR = $(BASEDIR)/test
DOCDIR= $(BASEDIR)/docs
OBJ_EXT=o
TEST_EXT=cpp

# Build a list of srcs and bins to build
SRCS=$(wildcard $(BASEDIR)/src/*.cc)
EXES=$(wildcard $(BASEDIR)/test/*.cpp)
OBJS=$(subst $(SRCDIR), $(OBJDIR),$(subst cc,$(OBJ_EXT),$(SRCS)))
BINS=$(subst $(TESTDIR), $(EXEDIR),$(subst .$(TEST_EXT),,$(EXES)))

.PHONY: all
all: lib $(BINS)

docs: all
	doxygen Doxyfile

$(EXEDIR)/%:  $(TESTDIR)/%.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/interface/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(OBJDIR)/%.$(OBJ_EXT):  $(SRCDIR)/%.cc $(BASEDIR)/interface/%.h*
	$(CXX) $(CXXFLAGS) -fPIC -c $<  -o $@

$(LIBDIR)/lib$(LIBNAME).so:  $(OBJS) $(LIBDEPENDS) 
	$(LD) $(LDFLAGS) -o $(LIBDIR)/lib$(LIBNAME).so $(OBJS) $(LIBS)

lib: $(LIBDIR)/lib$(LIBNAME).so

# info:
# 	@echo "LIBS: " $(LIBS)
# 	@echo "CXXFLAGS: " $(CXXFLAGS)
# 	@echo "Source files: " $(SRCS) 
# 	@echo "Object files: " $(OBJS)
# 	@echo "Executables:  " $(TARGETS)

clean:
	rm -rf $(OBJS) $(LIBDIR)/lib$(LIBNAME).so $(BINS)

