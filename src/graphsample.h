#include <string>

#ifndef GRAPHSAMPLE_H_
#define GRAPHSAMPLE_H_

class GraphSample {

  std::string left, right, output;
  int k;
  double rate;

public:

    GraphSample(
      std::string l,
      std::string r,
      std::string o,
      int k_,
      double rate_
    ) : left(l), right(r), output(o), k(k_), rate(rate_) { };

    void run(int usrseed, bool diginorm, bool only_part);
};

#endif  // GRAPHSAMPLE_H_
