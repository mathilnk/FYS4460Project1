#ifndef THERMOSTATTEST_H
#define THERMOSTATTEST_H
#include"cellsolver.h"
class ThermostatTest
{
public:
    ThermostatTest();
    void runBandA_Always();
    void runBUpDown(double first, double second, double first_timestep, double second_timestep);
};

#endif // THERMOSTATTEST_H
