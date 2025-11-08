#include "output_writer/write_paraview_output.h"

#include <iostream>
#include <cstdlib>
#include "settings/settings.h"

int main(int argc, char *argv[])
{
  // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
  if (argc == 1)
  {
    std::cout << "usage: " << argv[0] << " <filename>" << std::endl;

    return EXIT_FAILURE;
  }

  // read in the first argument
  std::string filename = argv[1];

  std::cout << "Parsing settings..." << std::endl;
#ifndef NDEBUG
  // print message
  std::cout << "filename: \"" << filename << "\"" << std::endl;
#endif

  SettingsParser settingsParser;
  Settings settings = settingsParser.loadFromFile(filename);
  std::cout << "Finished parsing!" << std::endl;

#ifndef NDEBUG
  std::cout << settings;
#endif

  return EXIT_SUCCESS;
}
