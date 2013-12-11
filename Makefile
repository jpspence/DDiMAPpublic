include ./findcudalib.mk

# Location of the CUDA Toolkit
CUDA_PATH       ?= /usr/local/cuda-5.5

# internal flags
NVCCFLAGS   := -m${OS_SIZE}
CCFLAGS     :=
NVCCLDFLAGS :=
LDFLAGS     :=

# Extra user flags
EXTRA_NVCCFLAGS   ?=
EXTRA_NVCCLDFLAGS ?=
EXTRA_LDFLAGS     ?= 
EXTRA_CCFLAGS     ?=

# OS-specific build flags
ifneq ($(DARWIN),) 
  LDFLAGS += -rpath $(CUDA_PATH)/lib
  CCFLAGS += -arch $(OS_ARCH) $(STDLIB)
else
  ifeq ($(OS_ARCH),armv7l)
    ifeq ($(abi),gnueabi)
      CCFLAGS += -mfloat-abi=softfp
    else
      # default to gnueabihf
      override abi := gnueabihf
      LDFLAGS += --dynamic-linker=/lib/ld-linux-armhf.so.3
      CCFLAGS += -mfloat-abi=hard
    endif
  endif
endif

ifeq ($(ARMv7),1)
NVCCFLAGS += -target-cpu-arch ARM
ifneq ($(TARGET_FS),) 
CCFLAGS += --sysroot=$(TARGET_FS)
LDFLAGS += --sysroot=$(TARGET_FS)
LDFLAGS += -rpath-link=$(TARGET_FS)/lib
LDFLAGS += -rpath-link=$(TARGET_FS)/usr/lib
LDFLAGS += -rpath-link=$(TARGET_FS)/usr/lib/arm-linux-$(abi)
endif
endif

# Debug build flags
ifeq ($(dbg),1)
      NVCCFLAGS += -g -G
      TARGET := debug
else
      TARGET := release
endif

ALL_CCFLAGS :=
ALL_CCFLAGS += $(NVCCFLAGS)
ALL_CCFLAGS += $(addprefix -Xcompiler ,$(CCFLAGS))
ALL_CCFLAGS += $(EXTRA_NVCCFLAGS)
ALL_CCFLAGS += $(addprefix -Xcompiler ,$(EXTRA_CCFLAGS))

ALL_LDFLAGS :=
ALL_LDFLAGS += $(ALL_CCFLAGS)
ALL_LDFLAGS += $(NVCCLDFLAGS)
ALL_LDFLAGS += $(addprefix -Xlinker ,$(LDFLAGS))
ALL_LDFLAGS += $(EXTRA_NVCCLDFLAGS)
ALL_LDFLAGS += $(addprefix -Xlinker ,$(EXTRA_LDFLAGS))

# Common includes and paths for CUDA
INCLUDES  := -Icommon/inc
LIBRARIES :=

################################################################################

# CUDA code generation flags
ifneq ($(OS_ARCH),armv7l)
GENCODE_SM10    := -gencode arch=compute_10,code=sm_10
endif
GENCODE_SM20    := -gencode arch=compute_20,code=sm_20
GENCODE_SM30    := -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=\"sm_35,compute_35\"
GENCODE_FLAGS   := $(GENCODE_SM10) $(GENCODE_SM20) $(GENCODE_SM30)

################################################################################

# Target rules


# Location of Bamtools on your machine. 
BAMTOOLS = include/bamtools

CXXFLAGS = -O2 -g -Wall -fmessage-length=0 -I$(BAMTOOLS)/src -I$(BAMTOOLS)
LIBS 	 = -L$(BAMTOOLS)/lib -lbamtools
OBJS 	 = bin/DDiMAP.o bin/DDiMAP-lib.o
TARGET   = DDiMAP

OBJST 	= bin/DDiMAP-threads.o bin/DDiMAP-lib.o
TARGETT = DDiMAP-threads 

all: $(TARGET) $(TARGETT) 

$(TARGETT): bin $(OBJST)
	$(CXX) -o bin/$(TARGETT) $(OBJST) $(LIBS)
	
$(TARGET): bin $(OBJS)
	$(CXX) -o bin/$(TARGET) $(OBJS) $(LIBS)

build: DDiMAPGPU

bin/DDiMAPGPU.o: src/DDiMAPGPU.cu
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $< -I$(BAMTOOLS)/src -I$(BAMTOOLS)

DDiMAPGPU: bin/DDiMAPGPU.o bin/DDiMAP-lib.o
	$(NVCC) $(ALL_LDFLAGS) -o bin/$@ $+ $(LIBRARIES) $(LIBS) 	

run: build
	./bin/DDiMAPGPU

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
