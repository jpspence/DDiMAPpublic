# Location of Bamtools on your machine. 
BAMTOOLS = include/bamtools
FASTQ    = include/readfq

CXXFLAGS = -O2 -g -Wall -fmessage-length=0 -I$(BAMTOOLS)/src -I$(FASTQ) -I$(BAMTOOLS) -std=c++0x
LIBS 	 = -L$(BAMTOOLS)/lib -lbamtools -lz

OBJS 	 = bin/DDiMAP.o bin/DDiMAP-lib.o bin/Read-Helper.o
TARGET   = DDiMAP

OBJST 	= bin/DDiMAP-threads.o bin/DDiMAP-lib.o bin/Read-Helper.o
TARGETT = DDiMAP-threads 

all: $(TARGET)

$(TARGET): bin $(OBJS)
	$(CXX) -o bin/$(TARGET) $(OBJS) $(LIBS)

bin:
	mkdir -p bin output

bin/%.o : src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

setup:
	git submodule update --init --recursive 
	cd include/bamtools && mkdir -p build && cd build && cmake .. && make && sudo make install && cd ../../..

clean:
	rm -rf bin/DDiMAP bin/*.o  output
