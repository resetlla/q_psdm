NVCC        := nvcc
CXX         := g++

SRC_DIR     := src
BIN_DIR     := bin
INC_DIR     := $(SRC_DIR)

TARGET      := $(BIN_DIR)/q-psdm

CU_SRCS     := $(SRC_DIR)/q-psdm.cu
CC_SRCS     := $(SRC_DIR)/fft.cc $(SRC_DIR)/segy.cc $(SRC_DIR)/psdmpkg.cc
OBJS        := $(CU_SRCS:.cu=.o) $(CC_SRCS:.cc=.o)

NVCCFLAGS   := -O3 -arch=sm_70 -I$(INC_DIR)
CXXFLAGS    := -O2 -I$(INC_DIR)
LDFLAGS     := 

all: $(TARGET)

$(TARGET): $(OBJS) | $(BIN_DIR)
	$(NVCC) $(NVCCFLAGS) $(OBJS) -o $@ $(LDFLAGS)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

run: $(TARGET)
	$(TARGET)

clean:
	rm -rf $(SRC_DIR)/*.o $(BIN_DIR)

.PHONY: all clean run