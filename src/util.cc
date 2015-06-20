#include "util.h"

string output_path_from_input(string &input, string &output) {
  int lastComponentIdx = input.rfind('/');
  string filename = input.substr(lastComponentIdx + 1);
  return output + "." + filename;
}

bool file_exists (const string& name) {
  ifstream f(name.c_str());
  if (f.good()) {
    // f.close();
    return true;
  } else {
    // f.close(); 
    return false;
  }
}
