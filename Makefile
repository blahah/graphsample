LIBKHMER_DIR = "./third-party/khmer/lib"

all: libkhmer
	mkdir -p bin
	g++ -std=c++11 -o bin/graphsample graphsample.cpp -I $(LIBKHMER_DIR) $(LIBKHMER_DIR)/libkhmer.a -I $(TCLAP_DIR)/include

libkhmer:
	$(MAKE) -C $(LIBKHMER_DIR)
