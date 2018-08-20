#====== for tfarm =======#
#CC=/usr/local/bin/gcc
#CXX=/usr/local/bin/g++
#ROOT	= /home/dragon/sfw/root-6.08.06-install/bin/root
#====== for farm =======#
#CC=gcc
#CXX=g++
#ROOT	= /usr/local/cern/root_v5.32.04.x86_64_fc8_gcc4.1.2/bin
#====== for gaio machine =======#
CC	= /usr/local/gcc-8.1.0/bin/gcc
C++	= /usr/local/gcc-8.1.0/bin/g++
ROOT	= /usr/local/hep/root_v6.14.00/bin/root

#CC=icc
#CXX=icpc
#CFLAGS  = -parallel -par-report3 -par-threshold0 -O3

CFLAGS  = -O2

BINDIR = ./bin
LIBDIR = ./lib

ROOTFLAGS = $(shell $(ROOT)-config --cflags)
ROOTLIBS = $(shell $(ROOT)-config --libs) -lMinuit
ROOTGLIBS = $(shell $(ROOT)-config --glibs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS)

TARGET1=     AC_ana1
OBJS1=       AC_ana1.o ParamMan.o Tree.o Setting.o
TARGET2=    Pedestal_ana
OBJS2=      Pedestal_ana.o Tree.o Setting.o
TARGET3=     AC_ana1_itabashi
OBJS3=       AC_ana1_itabashi.o ParamMan.o Tree.o

all: $(TARGET1) \
     $(TARGET2) \

$(LIBDIR)/%.o : %.cc
	$(CXX) $(CFLAGS) -c -o $@ $< $(CXXFLAGS)

$(TARGET1): $(patsubst %,$(LIBDIR)/%,$(OBJS1))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET2): $(patsubst %,$(LIBDIR)/%,$(OBJS2))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

$(TARGET3): $(patsubst %,$(LIBDIR)/%,$(OBJS2))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS)

.PHONY: clean
clean:
	rm -f $(LIBDIR)/*.o core $(BINDIR)/*

