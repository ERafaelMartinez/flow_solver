#pragma once

#include <iostream>
#include <array>

/** All settings that parametrize a simulation run.
 */
struct Settings
{
    std::array<int, 2> nCells;          //< number of cells in x and y direction
    std::array<double, 2> physicalSize; //< physical size of the domain

    double re = 1000;       //< reynolds number
    double endTime = 10.0;  //< end time of the simulation
    double tau = 0.5;       //< safety factor for time step width
    double maximumDt = 0.1; //< maximum time step width

    std::array<double, 2> g{0., 0.}; //< external forces

    bool useDonorCell = false; //< if the donor cell scheme schould be used
    double alpha = 0.5;        //< factor for donor-cell scheme

    std::array<double, 2> dirichletBcBottom; //< prescribed values of u,v at bottom of domain
    std::array<double, 2> dirichletBcTop;    //< prescribed values of u,v at top of domain
    std::array<double, 2> dirichletBcLeft;   //< prescribed values of u,v at left of domain
    std::array<double, 2> dirichletBcRight;  //< prescribed values of u,v at right of domain

    std::string pressureSolver = "SOR"; //< which pressure solver to use, "GaussSeidel" or "SOR"

    double omega = 1.0;                  //< overrelaxation factor
    double epsilon = 1e-5;               //< tolerance for the residual in the pressure solver
    int maximumNumberOfIterations = 1e5; //< maximum number of iterations in the solver

    /// @brief writes a string representation into the output stream
    friend std::ostream &operator<<(std::ostream &os, const Settings &s);
};

class SettingsParser
{
public:
    //! parse a text file with settings, each line contains "<parameterName> = <value>"
    Settings loadFromFile(std::string filename);

private:
    void parseParameter(Settings &settings, std::string parameter_name, std::string value);
};