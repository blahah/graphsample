#include <string>

#include "hashbits.hh"
#include "subset.hh"

#ifndef PARTITION_H_
#define PARTITION_H_

using namespace std;

namespace khmer {

  class Partition : public SubsetPartition {
  public:
    Partition(Hashtable * ht) : SubsetPartition(ht) {};

    size_t output_sampled_partitions(
      const string &left,
      const string &right,
      const string &out_left,
      const string &out_right,
      double rate,
      int usrseed
    );
  };

}

#endif  // PARTITION_H_
