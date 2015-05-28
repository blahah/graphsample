#include "graphsample.h"

#include <string>
#include <tclap/CmdLine.h>

using namespace TCLAP;

int main (int argc, char* argv[]) {

  try {

    // parse arguments

    CmdLine cmd("graphsample - representative FASTQ sampling",
                ' ', "0.0.1");

    ValueArg<std::string> leftArg("l", "left", "Left read file in FASTQ format", true, "", "string", cmd);

    ValueArg<std::string> rightArg("r", "right", "Right read file in FASTQ format", true, "", "string", cmd);

    ValueArg<std::string> outputArg("o", "output", "Output file prefix", false, "sampled", "string", cmd);

    ValueArg<int> kArg("k", "wordsize", "Word size for building de-Bruijn graph", false, 19, "int", cmd);

    ValueArg<double> rateArg("a", "rate", "Rate at which to sample components", false, 0.2, "double", cmd);

    ValueArg<int> seedArg("s", "seed", "Random seed >= 0 (default uses time)", false, -1, "int", cmd);

    cmd.parse(argc, argv);

    std::string left = leftArg.getValue();
    std::string right = rightArg.getValue();
    std::string output = outputArg.getValue();
    int k = kArg.getValue();
    double rate = rateArg.getValue();
    int usrseed = seedArg.getValue();

    // run program

    GraphSample gs(left, right, output, k, rate);
    gs.run(usrseed);

  } catch (ArgException &e) {
    std::cerr << "ERROR: " << e.error() << " for argument " <<
      e.argId() << std::endl;
  }

}
