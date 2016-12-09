#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H
#include <fstream>

class System; // Promise the compiler that this is a class even though we haven't included system.h here

class StatisticsSampler
{
private:
    std::ofstream m_file;
    double m_kineticEnergy = 0;
    double m_potentialEnergy = 0;
    double m_totalEnergy = 0;
    double m_temperature = 0;
    double m_density = 0;
    double m_diffusionConstant = 0;
    double m_volume = 0;
    double m_msd = 0;
public:
    StatisticsSampler();
    void saveToFile(System &system);
    void sample(System &system);
    void sampleKineticEnergy(System &system);
    void samplePotentialEnergy(System &system);
    void sampletotalEnergy(System &system);
    void sampleTemperature(System &system);
    void sampleDensity(System &system);
    void sampleMSD(System &system);
    void sampleDiffusionConstant(System &system);
    void samplevolume(System &system);

    double kineticEnergy() { return m_kineticEnergy; }
    double potentialEnergy() { return m_potentialEnergy; }
    double totalEnergy() { return m_kineticEnergy+m_potentialEnergy; }
    double temperature() { return m_temperature; }
    double density() { return m_density; }
    double volume() {return m_volume;}
    double msd() {return m_msd;}
    double diffusionConstant() {return m_diffusionConstant;}
};
#endif
