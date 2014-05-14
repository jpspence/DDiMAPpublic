# Common includes and paths for CUDA
INCLUDES  := -Icommon/inc
LIBRARIES :=

# Location of Bamtools on your machine. 
BAMTOOLS = include/bamtools
FASTQ    = include/readfq

CXXFLAGS = -O2 -g -Wall -fmessage-length=0 -I$(BAMTOOLS)/src -I$(FASTQ) -I$(BAMTOOLS) -std=c++0x
LIBS 	 = -L$(BAMTOOLS)/lib -lbamtools -lz

OBJS 	 = bin/DDiMAP.o bin/DDiMAP-lib.o bin/DDiMAP-test.o bin/Read-Helper.o # bin/gnuplot_i.o 
TARGET   = DDiMAP

OBJST 	= bin/DDiMAP-threads.o bin/DDiMAP-lib.o bin/Read-Helper.o
TARGETT = DDiMAP-threads 

all: $(TARGET) #  $(TARGETT)  DDiMAPGPU 

$(TARGETT): bin $(OBJST)
	$(CXX) -o bin/$(TARGETT) $(OBJST) $(LIBS)
	
$(TARGET): bin $(OBJS)
	$(CXX) -o bin/$(TARGET) $(OBJS) $(LIBS)

bin/gnuplot_i.o: include/gnuplot_i/gnuplot_i.c include/gnuplot_i/gnuplot_i.h
	gcc -c -o bin/gnuplot_i.o include/gnuplot_i/gnuplot_i.c

bin:
	mkdir -p bin

bin/%.o : src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

setup:
	git submodule update --init --recursive 
	cd include/bamtools && mkdir -p build && cd build && cmake .. && make && sudo make install && cd ../../..
	make -C include/cu

clean:
	rm -rf bin
