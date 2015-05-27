#include <string>

#include "hashbits.hh"
#include "subset.hh"

#ifndef PARTITION_H_
#define PARTITION_H_

namespace khmer {

  class Partition : public SubsetPartition {
  public:
    Partition(Hashtable * ht) : SubsetPartition(ht) {};
    size_t output_sampled_partitions(const std::string &infilename,
                              const std::string &outputfile, double rate);
  };

}

#endif  // PARTITION_H_
