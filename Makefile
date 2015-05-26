LIBKHMER_DIR = "./third-party/khmer/lib"

all: libkhmer
	mkdir -p bin
	g++ -std=c++11 -o bin/graphsample graphsample.cpp -I ./third-party/khmer/lib -L ./third-party/khmer/lib/ -lkhmer

libkhmer:
	$(MAKE) -C $(LIBKHMER_DIR)
