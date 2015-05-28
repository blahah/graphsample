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

void GraphSample::run(void) {

  // for now we set in stone, min_tables = 2 and hashsize = 4e9
  // we therefore want two primes > 4e9,
  // the smallest of which are 4e9+7 and 4e9+9
  vector<HashIntoType> sizes =
    {static_cast<HashIntoType>(4e9+7), static_cast<HashIntoType>(4e9+9)};

  Hashbits *htable = new Hashbits(k, sizes);

  unsigned long long n_consumed = 0;
  unsigned int total_reads = 0;
  htable->consume_fasta_and_tag(left, total_reads, n_consumed);
  htable->consume_fasta_and_tag(right, total_reads, n_consumed);

  // partition
  Partition part(htable);
  part.do_partition(0, 0, false, true);

  // subsample partitions and output reads
  string out_left = output_path_from_input(left, output);
  string out_right = output_path_from_input(right, output);
  part.output_sampled_partitions(left, right,
                                 out_left, out_right, rate);

}
