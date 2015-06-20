#ifndef UTIL_H_
#define UTIL_H_

#include <string>
#include <fstream>

using namespace std;

extern string output_path_from_input(string &input, string &output);

extern bool file_exists (const string& name);

#endif  // UTIL_H_
