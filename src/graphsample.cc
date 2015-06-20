#include "graphsample.h"
#include "util.h"

#include <vector>
#include <random>
#include <algorithm>
#include <exception>

#include "partition.h"
#include "hashbits.hh"

using namespace khmer;
using namespace std;

void GraphSample::run(int usrseed, bool diginorm, bool only_part) {

  cout << "Building kmer hash table from reads" << endl;

  // for now we set in stone, min_tables = 2 and hashsize = 4e9
  // we therefore want two primes > 4e9,
  // the smallest of which are 4e9+7 and 4e9+9
  vector<HashIntoType> sizes =
    {static_cast<HashIntoType>(4e9+7), static_cast<HashIntoType>(4e9+9)};

  Hashbits htable(k, sizes);
  Partition part(&htable);

  string htableprefix = "countinghash";
  string htable_out = output_path_from_input(left, htableprefix) + ".htable";

  if (file_exists(htable_out)) {

    cout << "Loading counting hash map from existing file" << endl;

    try {
      htable.load(htable_out);
    } catch (exception &e) {
      cerr << "ERROR loading counting hash map" << endl << e.what() << endl;
      exit(1);
    }

  } else {

    // load reads
    unsigned long long n_consumed = 0;
    unsigned int total_reads = 0;
    htable.consume_fasta_and_tag(left, total_reads, n_consumed);
    htable.consume_fasta_and_tag(right, total_reads, n_consumed);
    cout << "Consumed " << n_consumed << " kmers from " << total_reads
         << " read pairs" << endl;

    htable.save(htable_out);

  }


  string pmapprefix = "partitioned";
  string pmap_out = output_path_from_input(left, pmapprefix) + ".pmap";

  if (file_exists(pmap_out)) {

    // load existing pmap

    cout << "Loading partition map from existing file" << endl;

    try {
      part.load_partitionmap(pmap_out);
    } catch (exception &e) {
      cerr << "ERROR loading partition map" << endl << e.what() << endl;
      exit(1);
    }

    cout << "Partitions loaded" << endl;

    size_t n_partitions = 0;
    size_t n_unassigned = 0;
    part.count_partitions(n_partitions, n_unassigned);
    cout << "Loaded " << n_partitions << " partitions from saved map" << endl;

  } else  {

    // partition
    cout << "Partitioning graph" << endl;
    part.do_partition(0, 0, true, true);
    size_t n_partitions = 0;
    size_t n_unassigned = 0;
    part.count_partitions(n_partitions, n_unassigned);
    cout << "Found " << n_partitions << " partitions" << endl;

    // // join by read pairing
    // cout << "Joining partitions" << endl;
    // int n_joined = part.join_bridged_partitions(left, right);
    // cout << "Joined " << n_joined << " partitions" << endl;

    // save
    part.save_partitionmap(pmap_out);

  }

  if (only_part) {

    string prefix = "partitioned";
    string out_left = output_path_from_input(left, prefix);
    string out_right = output_path_from_input(right, prefix);
    cout << out_left;
    if (!file_exists(out_left) && !file_exists(out_right)) {
      part.output_partitions(left, right, out_left, out_right);
    } else {
      cout << "Partitioned read files already exist! Not overwriting" << endl;
    }

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
