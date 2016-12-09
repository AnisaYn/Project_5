#include "system.h"
#include "statisticssampler.h"
#include "lennardjones.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>
using namespace std;

using std::ofstream; using std::cout; using std::endl;

StatisticsSampler::StatisticsSampler()
{

}

void StatisticsSampler::saveToFile(System &system)
{
    // Save the statistical properties for each timestep for plotting etc.
    // First, open the file if it's not open already
    if(!m_file.good()) {
        m_file.open("statistics.txt", ofstream::out);
        // If it's still not open, something bad happened...
        if(!m_file.good()) {
            cout << "Error, could not open statistics.txt" << endl;
            exit(1);
        }
        m_file << setw(20) << "Timestep" <<
                  setw(20) << "Time" <<
                  setw(20) << "KineticEnergy" <<
                  setw(20) << "PotentialEnergy" <<
                  setw(20) << "TotalEnergy" <<
                  setw(20) << "Temperature" <<
                  setw(20) << "Density" <<
                  setw(20) << "MSD" <<
                  setw(20) << "Diffusion Constant" <<
                  setw(20) << endl;

    }


    //  m_file << system.steps() << "  " << system.time() << "  " << m_kineticEnergy << "  " <<m_potentialEnergy<< "  " <<m_totalEnergy<< "  " <<m_temperature<< "  " <<m_density<< "  " <<m_diffusionConstant<< "  " << endl;
    // Print out values here

    m_file <<  system.steps() << " "
           << UnitConverter::timeToSI(system.time()) << " "
           << UnitConverter::energyToEv(kineticEnergy()) << " "
           << UnitConverter::energyToEv(potentialEnergy()) << " "
           << UnitConverter::energyToEv(totalEnergy()) << " "
           << UnitConverter::temperatureToSI(temperature()) << " "
           << density() << " "
           << msd() << " "
           << UnitConverter::diffusionToSI(diffusionConstant()) << " "// m^2/s
              << endl;
/*
    m_file << setw(20) <<  system.steps() <<
              setw(20) << UnitConverter::timeToSI(system.time()) <<
              setw(20) << UnitConverter::energyToEv(kineticEnergy()) <<
              setw(20) << UnitConverter::energyToEv(potentialEnergy()) <<
              setw(20) << UnitConverter::energyToEv(totalEnergy()) <<
              setw(20) << UnitConverter::temperatureToSI(temperature()) <<
              setw(20) << density() <<
              setw(20) << msd() <<
              setw(20) << UnitConverter::diffusionToSI(diffusionConstant()) << // m^2/s
              setw(20) << endl;
*/
}


void StatisticsSampler::sample(System &system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    sampleMSD(system);
    sampleDiffusionConstant(system);
    saveToFile(system);
}

void StatisticsSampler::sampleKineticEnergy(System &system)
{
    m_kineticEnergy = 0; // Remember to reset the value from the previous timestep
    for(Atom *atom : system.atoms()) {
        m_kineticEnergy += 0.5*atom->mass()*atom->velocity.lengthSquared();
    }

}

void StatisticsSampler::samplePotentialEnergy(System &system)
{
    m_potentialEnergy = system.potential().potentialEnergy();

}

void StatisticsSampler:: sampletotalEnergy(System &system)
{
    m_totalEnergy = m_kineticEnergy+m_potentialEnergy;
}

void StatisticsSampler::sampleTemperature(System &system)
{
    m_temperature = (2.0/3.0)*m_kineticEnergy/system.atoms().size();// Hint: reuse the kinetic energy that we already calculated

}

//void StatisticsSampler::samplevolume(System &system)
//{
//  m_volume = atom->mass()/system.atoms().size();
//}

void StatisticsSampler::sampleDensity(System &system)
{

    m_density = system.atoms().size()/system.volume();

}

void StatisticsSampler::sampleMSD(System &system)
{
    m_msd = 0;
    for(Atom *atom : system.atoms()){
        m_msd += (atom->position - atom->initial_position).lengthSquared();

    }
    m_msd /= system.atoms().size();
}

void StatisticsSampler::sampleDiffusionConstant(System &system)
{
    m_diffusionConstant = m_msd/6/system.time();

}
