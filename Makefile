# Location of Bamtools on your machine. 
BAMTOOLS = include/bamtools

CXXFLAGS =	-O2 -g -Wall -fmessage-length=0 -I$(BAMTOOLS)/src -I$(BAMTOOLS)
LIBS 		= -L$(BAMTOOLS)/lib -lbamtools
OBJS 		=	bin/DDiMAP.o bin/DDiMAP-lib.o
TARGET  =	DDiMAP

OBJST 	= bin/DDiMAP-threads.o bin/DDiMAP-lib.o
TARGETT = DDiMAP-threads 

all: $(TARGET) $(TARGETT)

$(TARGETT): bin $(OBJST)
	$(CXX) -o bin/$(TARGETT) $(OBJST) $(LIBS)
	
$(TARGET): bin $(OBJS)
	$(CXX) -o bin/$(TARGET) $(OBJS) $(LIBS)


bin:
	mkdir -p bin

bin/%.o : src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

install:
	git submodule update --init --recursive 
	cd include/bamtools && mkdir -p build && cd build && cmake .. && make && cd ../../..
	make -C include/cu

clean:
	rm -rf bin
