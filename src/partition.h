#include <string>

#include "hashbits.hh"
#include "subset.hh"
#include "counting.hh"
#include "read_parsers.hh"

#ifndef PARTITION_H_
#define PARTITION_H_

using namespace std;

namespace khmer {

  class Partition : public SubsetPartition {
  public:
    Partition(Hashtable * ht) : SubsetPartition(ht) {};

    size_t output_partitions(
      const string &left,
      const string &right,
      const string &out_left,
      const string &out_right
    );

    size_t output_sampled_partitions(
      const string &left,
      const string &right,
      const string &out_left,
      const string &out_right,
      double rate,
      int usrseed,
      int k,
      bool diginorm
    );

    size_t join_bridged_partitions(
        const string &left,
        const string &right
    );

    bool below_min_coverage(
      read_parsers::Read& read,
      CountingHash& hash,
      int k
    );

  };

}

#endif  // PARTITION_H_
