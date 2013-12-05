# Location of Bamtools on your machine. 
BAMTOOLS = bamtools

CXXFLAGS =	-O2 -g -Wall -fmessage-length=0 -I$(BAMTOOLS)/src -I$(BAMTOOLS)
LIBS = -L$(BAMTOOLS)/lib -lbamtools
OBJS =		DDiMAP.o
TARGET =	DDiMAP

	
$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all: $(TARGET)
	
install:
	git clone git@github.com:pezmaster31/bamtools.git
	cd bamtools && mkdir -p build && cd build && cmake .. && make

clean:
	rm -f $(OBJS) $(TARGET)
