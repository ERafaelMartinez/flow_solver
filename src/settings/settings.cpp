#include <fstream>   // for file operations
#include <iostream>  // for cout

void loadFromFile(std::string filename)
{
    // open the file
    std::ifstream file(filename.c_str(), std::ios::in);

    // check if file is open
    if (!file.is_open())
    {
        std::cout << "Error: Could not open parameter file \"" << filename << "\"." << std::endl;
        return;
    }


    // loop over each line in the file
    for (int lineNo = 0; !file.eof(); lineNo++)
    {
        std::string line;
        std::getline(file, line);

        // print the line
        std::cout << "Line " << lineNo << ": " << line << std::endl;
    }
}