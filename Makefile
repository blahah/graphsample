LIBKHMER_DIR = "./third-party/khmer/lib"
TCLAP_DIR = "./third-party/tclap-1.2.1"

all: libkhmer
	mkdir -p bin
	g++ -std=c++11 -o bin/graphsample src/*.cc -I $(LIBKHMER_DIR) $(LIBKHMER_DIR)/libkhmer.a -I $(TCLAP_DIR)/include

debug: libkhmer
	mkdir -p bin
	g++ -std=c++11 -g -o bin/graphsample_dbg src/*.cc -I $(LIBKHMER_DIR) $(LIBKHMER_DIR)/libkhmer.a -I $(TCLAP_DIR)/include

libkhmer:
	$(MAKE) -C $(LIBKHMER_DIR)
