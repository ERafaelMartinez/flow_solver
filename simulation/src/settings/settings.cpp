#include "settings.h"

#include <iostream>

std::ostream &operator<<(std::ostream &os, const Settings &s)
{
    os << "Settings: " << std::endl
       << "  physicalSize: " << s.physicalSize[0] << " x " << s.physicalSize[1] << ", nCells: " << s.nCells[0] << " x " << s.nCells[1] << std::endl
       << "  endTime: " << s.endTime << " s, re: " << s.re << ", g: (" << s.g[0] << "," << s.g[1] << "), tau: " << s.tau << ", maximum dt: " << s.maximumDt << std::endl
       << "  dirichletBC: " << std::endl
       << "    bottom: (" << s.dirichletBcBottom[0] << "," << s.dirichletBcBottom[1] << ")" << std::endl
       << "    top: (" << s.dirichletBcTop[0] << "," << s.dirichletBcTop[1] << ")" << std::endl
       << "    left: (" << s.dirichletBcLeft[0] << "," << s.dirichletBcLeft[1] << ")" << std::endl
       << "    right: (" << s.dirichletBcRight[0] << "," << s.dirichletBcRight[1] << ")" << std::endl
       << "  useDonorCell: " << std::boolalpha << s.useDonorCell << std::endl
       << "    alpha: " << s.alpha << std::endl
       << "  pressureSolver: " << s.pressureSolver << std::endl
       << "    omega: " << s.omega << std::endl
       << "    epsilon: " << s.epsilon << std::endl
       << "    maximumNumberOfIterations: " << s.maximumNumberOfIterations << std::endl;
    return os;
}