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

void GraphSample::run(int usrseed, bool diginorm, bool only_part) {

  cout << "Building kmer hash table from reads" << endl;
  // for now we set in stone, min_tables = 2 and hashsize = 4e9
  // we therefore want two primes > 4e9,
  // the smallest of which are 4e9+7 and 4e9+9
  vector<HashIntoType> sizes =
    {static_cast<HashIntoType>(4e9+7), static_cast<HashIntoType>(4e9+9)};

  Hashbits htable(k, sizes);

  unsigned long long n_consumed = 0;
  unsigned int total_reads = 0;
  htable.consume_fasta_and_tag(left, total_reads, n_consumed);
  htable.consume_fasta_and_tag(right, total_reads, n_consumed);
  cout << "Consumed " << n_consumed << " kmers from " << total_reads
       << " read pairs" << endl;

  // partition
  Partition part(&htable);
  cout << "Partitioning graph" << endl;
  part.do_partition(0, 0, true, true);
  size_t n_partitions = 0;
  size_t n_unassigned = 0;
  part.count_partitions(n_partitions, n_unassigned);
  cout << "Found " << n_partitions << " partitions" << endl;

  // join by read pairing
  cout << "Joining partitions" << endl;
  int n_joined = part.join_bridged_partitions(left, right);
  cout << "Joined " << n_joined << " partitions" << endl;

  if (only_part) {
    string prefix = "partitioned";
    string out_left = output_path_from_input(left, prefix);
    string out_right = output_path_from_input(right, prefix);
    part.output_partitions(left, right, out_left, out_right);

  } else {

    // subsample partitions and output reads
    cout << "Subsampling partitions and reads" << endl;
    string out_left = output_path_from_input(left, output);
    string out_right = output_path_from_input(right, output);
    part.output_sampled_partitions(left, right, out_left, out_right,
                                   rate, usrseed, k, diginorm);

  }

  cout << "Done!" << endl;
}
