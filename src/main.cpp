#include "settings/settings.h"
#include "simulation.h"
#include <cstdlib>
#include <iostream>
#include <mpi.h>

int main(int argc, char *argv[]) {
  // if the number of given command line arguments is only 1 (= the program
  // name), print out usage information and exit
  if (argc == 1) {
    std::cout << "usage: " << argv[0] << " <filename>" << std::endl;

    return EXIT_FAILURE;
  }

  MPI_Init(&argc, &argv);

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

  Simulation simuation(&settings);
  simuation.run();

  MPI_Finalize();
  return EXIT_SUCCESS;
}
