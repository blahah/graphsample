LIBKHMER_DIR = "./third-party/khmer/lib"

all: libkhmer
	g++ -std=c++11 -o graphsample graphsample.cpp -I ./third-party/lib -L ./third-party/lib/ -lkhmer

libkhmer:
	$(MAKE) -C $(LIBKHMER_DIR)
