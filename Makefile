# Location of Bamtools on your machine. 
BAMTOOLS = include/bamtools

CXXFLAGS =	-O2 -g -Wall -fmessage-length=0 -I$(BAMTOOLS)/src -I$(BAMTOOLS)
LIBS = -L$(BAMTOOLS)/lib -lbamtools
OBJS =		DDiMAP.o
TARGET =	DDiMAP

	
$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all: install $(TARGET)
	
install:
	git submodule update --init --recursive 
	cd include/bamtools && mkdir -p build && cd build && cmake .. && make && cd ../../..
	make -C include/cu

clean:
	rm -f $(OBJS) $(TARGET)
