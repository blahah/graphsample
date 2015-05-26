#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <tclap/CmdLine.h>

#include "sample.cpp"
#include "hashbits.hh"

using namespace TCLAP;
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

  try {

    // parse arguments

    CmdLine cmd("graphsample - representative FASTQ subsampling",
                ' ', "0.0.1");

    ValueArg<std::string> inputArg("i", "input", "Read file in interleaved FASTQ format", true, "", "string");
    cmd.add(inputArg);

    ValueArg<std::string> outputArg("o", "output", "Path to write subsampled reads to", true, "", "string");
    cmd.add(outputArg);

    cmd.parse(argc, argv);

    std::string input = inputArg.getValue();
    std::string output = outputArg.getValue();

    // run program

    Hashbits *htable = new_hashbits();

    // load reads into the table
    std::string read_pairs_in(input);
    unsigned long long n_consumed = 0;
    unsigned int total_reads = 0;
    htable->consume_fasta_and_tag(read_pairs_in, total_reads, n_consumed);

    // partition
    Part *subset_p = NULL;
    subset_p = new Part(htable);
    subset_p->do_partition(0, 0, false, true);

    // subsample partitions and output reads
    std::string read_pairs_out(output);
    subset_p->output_sampled_partitions(read_pairs_in, read_pairs_out);


  } catch (ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " <<
      e.argId() << std::endl;
  }

}
