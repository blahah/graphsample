// load a pair of fastq files into a compressible graph format
// partition the graph
// take a random subsample of the partitions
// align reads to the sampled partitions and keep those that align

#include <string>
#include <vector>
#include <random>

#include "sample.cpp"
#include "hashbits.hh"

using namespace khmer;

Hashbits* new_hashbits() {

  // for now we set in stone, min_tables = 2 and hashsize = 4e9
  // we therefore want two primes > 4e9,
  // the smallest of which are 4e9+7 and 4e9+9
  std::vector<HashIntoType> sizes =
    {static_cast<HashIntoType>(4e9+7), static_cast<HashIntoType>(4e9+9)};

  // kmer size
  int k = 31;
  Hashbits *htable = new Hashbits(k, sizes);

  return htable;

}

int main (int argc, char* argv[]) {

  Hashbits *htable = new_hashbits();

  // load reads into the table
  std::string read_pairs_in("interleaved.fq");
  unsigned long long n_consumed = 0;
  unsigned int total_reads = 0;
  htable->consume_fasta_and_tag(read_pairs_in, total_reads, n_consumed);

  // partition
  Part *subset_p = NULL;
  subset_p = new Part(htable);
  subset_p->do_partition(0, 0, false, true);

  // subsample partitions and output reads
  std::string read_pairs_out("partitioned.fq");
  subset_p->output_sampled_partitions(read_pairs_in, read_pairs_out);

}
