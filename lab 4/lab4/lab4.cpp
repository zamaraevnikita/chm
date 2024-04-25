#include <iostream>
#include "NonlinearEquations.h"

using namespace luMath;

int main()
{
    NonlinearEquations<double> data;
    switch (data.getMethod())
    {
    case NonlinearEquations<double>::METHOD::NEWTON:  // Метод Ньютона (+ Метод Итераций)
        data.Newton();
        break;
    case NonlinearEquations<double>::METHOD::STEEPESTDESCENT: // Метод наискорейшего спуска
        data.SteepestDescent();
        break;
    }
    return 0;
}