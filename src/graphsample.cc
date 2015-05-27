#include "graphsample.h"

#include <vector>
#include <random>
#include <algorithm>

#include "partition.h"
#include "hashbits.hh"

using namespace khmer;
using namespace std;

string output_path_from_input(string &input, string &output) {
  int lastComponentIdx = input.rfind('/');
  string filename = input.substr(lastComponentIdx + 1);
  return output + "." + filename;
}

// Load reads into the counting hash
void load_reads(string reads_in, Hashbits *htable) {
  unsigned long long n_consumed = 0;
  unsigned int total_reads = 0;
  htable->consume_fasta_and_tag(reads_in, total_reads, n_consumed);
}

// Sample from a read file
void sample_reads(string reads_in, string out, Partition &part, double rate) {
  string reads_out = output_path_from_input(reads_in, out);
  part.output_sampled_partitions(reads_in, reads_out, rate);
}


void GraphSample::run(void) {

  // for now we set in stone, min_tables = 2 and hashsize = 4e9
  // we therefore want two primes > 4e9,
  // the smallest of which are 4e9+7 and 4e9+9
  vector<HashIntoType> sizes =
    {static_cast<HashIntoType>(4e9+7), static_cast<HashIntoType>(4e9+9)};

  Hashbits *htable = new Hashbits(k, sizes);

  load_reads(left, htable);
  if (right.length() > 0) {
    load_reads(right, htable);
  }

  // partition
  Partition part(htable);
  part.do_partition(0, 0, false, true);

  // subsample partitions and output reads
  sample_reads(left, output, part, rate);
  if (right.length() > 0) {
    sample_reads(right, output, part, rate);
  }

}
