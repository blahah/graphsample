#include "graphsample.h"

#include <string>
#include <tclap/CmdLine.h>

using namespace TCLAP;
using namespace std;

typedef unsigned int PartitionID;

class GraphSampleOutput : public StdOutput
{
  public:

    virtual void usage(CmdLineInterface& c)
    {
      cout << endl << "graphsample - representative FASTQ sampling:" << endl;
      cout << endl << "http://github.com/blahah/graphsample" << endl << endl;

      list<Arg*> args = c.getArgList();
      for (auto arg : args) {
        if (arg->getName() == "ignore_rest") {
          continue;
        }
        cout << "   " << arg->longID() << endl;
        cout << "      " << arg->getDescription() << endl;
        cout << endl;
      }
    }

    virtual void version(CmdLineInterface& c)
    {
      cout << "graphsample v0.0.1" << endl;
    }
};

int main (int argc, char* argv[]) {

  try {

    // parse arguments

    CmdLine cmd("graphsample - representative FASTQ sampling",
                ' ', "0.0.1");
    GraphSampleOutput helpmsg;
		cmd.setOutput( &helpmsg );

    ValueArg<int> seedArg("s", "seed", "Random seed >= 0 (default uses time)", false, -1, "int", cmd);

    ValueArg<double> rateArg("a", "rate", "Rate at which to sample components [0.0, 1.0] (default 0.2)  ", false, 0.2, "double", cmd);

    ValueArg<int> kArg("k", "wordsize", "Word size for building de-Bruijn graph (default 19)", false, 19, "int", cmd);

    ValueArg<string> outputArg("o", "output", "Output file prefix (default sample)", false, "sampled", "string", cmd);

    ValueArg<string> rightArg("r", "right", "Right read file in FASTQ format", true, "", "string", cmd);

    ValueArg<string> leftArg("l", "left", "Left read file in FASTQ format", true, "", "string", cmd);

    SwitchArg normArg("d", "diginorm", "Digitally normalise the reads in each cluster to a minimum coverage of 20", cmd);

    SwitchArg partArg("p", "only-part", "Skip sampling step, and just output reads tagged with their partition", cmd);


    cmd.parse(argc, argv);

    string left = leftArg.getValue();
    string right = rightArg.getValue();
    string output = outputArg.getValue();
    int k = kArg.getValue();
    double rate = rateArg.getValue();
    int usrseed = seedArg.getValue();
    bool diginorm = normArg.getValue();
    bool onlypart = partArg.getValue();

    // run program

    GraphSample gs(left, right, output, k, rate);
    gs.run(usrseed, diginorm, onlypart);

  } catch (ArgException &e) {
    cerr << "ERROR: " << e.error() << " for argument " <<
      e.argId() << endl;
  }

}
