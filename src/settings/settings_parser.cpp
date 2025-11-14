#include "settings.h"
#include <fstream>  // for file operations
#include <iostream> // for cout

// TODO: Move to its own util file
const char *whitespaceCharacters = " \t\r\n\f\v";
const char commentChar = '#';
const char keyValueDelimiter = '=';

std::string trim(std::string &str)
{
    // trim all whitespace characters on the left
    str.erase(0, str.find_first_not_of(whitespaceCharacters));
    // trim all whitespace character on the right
    str.erase(str.find_last_not_of(whitespaceCharacters) + 1);
    return str;
}

std::string removeCommentFromLine(std::string &line)
{
    size_t comment_index = line.find(commentChar);
    if (comment_index != std::string::npos)
    {
        line.erase(comment_index);
    }
    return line;
}

bool isStringWhitespaceOrEmpty(std::string &str)
{
    return str.find_first_not_of(whitespaceCharacters) == std::string::npos;
}

void SettingsParser::parseParameter(Settings &settings, std::string parameter_name, std::string value_to_parse)
{
    if (parameter_name == "nCellsX")
    {
        settings.nCells[0] = atoi(value_to_parse.c_str());
    }
    else if (parameter_name == "nCellsY")
    {
        settings.nCells[1] = atoi(value_to_parse.c_str());
    }
    else if (parameter_name == "physicalSizeX")
    {
        settings.physicalSize[0] = atof(value_to_parse.c_str());
    }
    else if (parameter_name == "physicalSizeY")
    {
        settings.physicalSize[1] = atof(value_to_parse.c_str());
    }
    else if (parameter_name == "re")
    {
        settings.re = atof(value_to_parse.c_str());
    }
    else if (parameter_name == "endTime")
    {
        settings.endTime = atof(value_to_parse.c_str());
    }
    else if (parameter_name == "tau")
    {
        settings.tau = atof(value_to_parse.c_str());
    }
    else if (parameter_name == "maximumDt")
    {
        settings.maximumDt = atof(value_to_parse.c_str());
    }
    else if (parameter_name == "gX")
    {
        settings.g[0] = atof(value_to_parse.c_str());
    }
    else if (parameter_name == "gY")
    {
        settings.g[1] = atof(value_to_parse.c_str());
    }
    else if (parameter_name == "useDonorCell")
    {
        // TODO: Check if this is robust
        // ane maybe use isstreamstring with boolalpha
        settings.useDonorCell = value_to_parse == "true";
    }
    else if (parameter_name == "alpha")
    {
        settings.alpha = atof(value_to_parse.c_str());
    }

    else if (parameter_name == "dirichletBottomX")
    {
        settings.dirichletBcBottom[0] = atof(value_to_parse.c_str());
    }
    else if (parameter_name == "dirichletBottomY")
    {
        settings.dirichletBcBottom[1] = atof(value_to_parse.c_str());
    }

    else if (parameter_name == "dirichletTopX")
    {
        settings.dirichletBcTop[0] = atof(value_to_parse.c_str());
    }
    else if (parameter_name == "dirichletTopY")
    {
        settings.dirichletBcTop[1] = atof(value_to_parse.c_str());
    }

    else if (parameter_name == "dirichletLeftX")
    {
        settings.dirichletBcLeft[0] = atof(value_to_parse.c_str());
    }
    else if (parameter_name == "dirichletLeftY")
    {
        settings.dirichletBcLeft[1] = atof(value_to_parse.c_str());
    }

    else if (parameter_name == "dirichletRightX")
    {
        settings.dirichletBcRight[0] = atof(value_to_parse.c_str());
    }
    else if (parameter_name == "dirichletRightY")
    {
        settings.dirichletBcRight[1] = atof(value_to_parse.c_str());
    }

    else if (parameter_name == "pressureSolver")
    {
        // TODO: Use enum insteam
        settings.pressureSolver = value_to_parse;
    }
    else if (parameter_name == "omega")
    {
        settings.omega = atof(value_to_parse.c_str());
    }
    else if (parameter_name == "epsilon")
    {
        settings.epsilon = atof(value_to_parse.c_str());
    }
    else if (parameter_name == "maximumNumberOfIterations")
    {
        settings.maximumNumberOfIterations = atof(value_to_parse.c_str());
    }
}

Settings SettingsParser::loadFromFile(std::string filename)
{
    // open the file
    std::ifstream file(filename.c_str(), std::ios::in);

    Settings settings;

    // check if file is open
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open parameter file \"" << filename << "\"." << std::endl;
        // TODO: Do we want to return a bool or throw an exception here so that the app aborts?
        throw std::invalid_argument("file could not be opened!");
    }

    // loop over each line in the file
    for (int lineNo = 0; !file.eof(); lineNo++)
    {
        std::string line;
        std::getline(file, line);

#ifndef NDEBUG
        std::cout << "Parsing line " << lineNo << ": " << line << std::endl;
#endif

        line = removeCommentFromLine(line);
        line = trim(line);
        if (line.empty())
        {
#ifndef NDEBUG
            std::cout << "\t empty line or comment, skipping line" << std::endl;
#endif
            continue;
        }

        // split name and value
        size_t delimiter_index = line.find(keyValueDelimiter);
        if (delimiter_index == std::string::npos)
        {
#ifndef NDEBUG
            std::cout << "\t no equal sign found, skipping line" << std::endl;
#endif
            continue;
        }

        std::string parameter_name = line.substr(0, delimiter_index);
        parameter_name = trim(parameter_name);

        std::string value_to_parse = line.substr(delimiter_index + 1);
        value_to_parse = trim(value_to_parse);

        if (isStringWhitespaceOrEmpty(parameter_name) || isStringWhitespaceOrEmpty(value_to_parse))
        {
#ifndef NDEBUG
            std::cout << "\t parameter name or value are empty, skipping line" << std::endl;
#endif
            continue;
        }

#ifndef NDEBUG
        std::cout << "\t parsing parameter: " << parameter_name << ", with value: " << value_to_parse << std::endl;
#endif
        parseParameter(settings, parameter_name, value_to_parse);
    }

    return settings;
}