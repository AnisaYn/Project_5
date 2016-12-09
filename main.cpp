#include "math/random.h"
#include "lennardjones.h"
#include "velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>

using namespace std;

int main(int numberOfArguments, char **argumentList)
{
    //Random::seed(3857358758);

    int numberOfUnitCells = 5;
    double initialTemperature = UnitConverter::temperatureFromSI(450); // measured in Kelvin
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms
    double dt = UnitConverter::timeFromSI(1e-14); // Measured in seconds.

    // If a first argument is provided, it is the number of unit cells
    if(numberOfArguments > 1) numberOfUnitCells = atoi(argumentList[1]);
    // If a second argument is provided, it is the initial temperature (measured in kelvin)
    if(numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));
    // If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
    if(numberOfArguments > 3) latticeConstant = UnitConverter::lengthFromAngstroms(atof(argumentList[3]));
    if(numberOfArguments > 4) dt = UnitConverter::timeFromSI(atof(argumentList[4]));


    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    // vector<double> temperatures = {100, 200, 300, 400, 500, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 700, 800, 900};
    /*
    for (double temperatures = UnitConverter::temperatureFromSI(100.0); temperatures < UnitConverter::temperatureFromSI(1000.0); temperatures += UnitConverter::temperatureFromSI(100.0)){
      System system;
      system.createFCCLattice(numberOfUnitCells, latticeConstant, temperatures);
    */
    System system;
    system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature);
    system.potential().setEpsilon(1.0);
    system.potential().setSigma(3.405);
    system.resetForce();
    system.removeTotalMomentum();

    StatisticsSampler statisticsSampler;

    IO movie("movie.xyz"); // To write the state to file

    cout << setw(20) << "Timestep" <<
            setw(20) << "Time" <<
            setw(20) << "Temperature" <<
            setw(20) << "KineticEnergy" <<
            setw(20) << "PotentialEnergy" <<
            setw(20) << "TotalEnergy" << endl;
    for(int timestep=0; timestep<50000; timestep++) {
        system.step(dt);
        statisticsSampler.sample(system);
        if( !(timestep % 100) ) {
            // Print the timestep every 100 timesteps
            cout << setw(20) <<  system.steps() <<
                    setw(20) << UnitConverter::timeToSI(system.time()) <<
                    setw(20) << UnitConverter::temperatureToSI(statisticsSampler.temperature()) <<
                    setw(20) << UnitConverter::energyToEv(statisticsSampler.kineticEnergy()) <<
                    setw(20) << UnitConverter::energyFromEv(statisticsSampler.potentialEnergy()) <<
                    setw(20) << UnitConverter::energyFromEv(statisticsSampler.totalEnergy()) <<
                    setw(20) << endl;

        }
        // movie.saveState(system);
    }

    movie.close();

    return 0;
}
