#include <iostream>
#include "Eigenvalue.h"
using namespace luMath;
int main()
{
    setlocale(LC_ALL, "Rus");
    Eigenvalue<double> matrix;

    switch (matrix.getTask())
    {
    case Eigenvalue<double>::TASK::EIGENVALUES:
        matrix.getEigenvalues();
        break;
    case Eigenvalue<double>::TASK::EIGENVECTORS:
        matrix.getEigenvalues();
        matrix.getEigenvectors();
        break;
    case Eigenvalue<double>::TASK::ROOTS:
        matrix.GetRoots();
        break;
    }

    return 0;
}

