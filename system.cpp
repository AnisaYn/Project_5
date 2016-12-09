#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"
#include <iomanip>
#include <iostream>

using std::cout;
using std::endl;
using std::setprecision;

System::System()
{

}

System::~System()
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
    for(Atom *atom : m_atoms){
        if(atom->position.x() > m_systemSize.x()){
            atom->position[0] -= m_systemSize[0];
            atom->initial_position[0] -= m_systemSize[0];
        }
        if(atom->position.y() > m_systemSize.y()){
            atom->position[1] -= m_systemSize[1];
            atom->initial_position[1] -= m_systemSize[1];
        }
        if(atom->position.z() > m_systemSize.z()){
            atom->position[2] -= m_systemSize[2];
            atom->initial_position[2] -= m_systemSize[2];
        }
        if(atom->position.x() < 0){
            atom->position[0] += m_systemSize[0];
            atom->initial_position[0] += m_systemSize[0];
        }
        if(atom->position.y() < 0){
            atom->position[1] += m_systemSize[1];
            atom->initial_position[1] += m_systemSize[1];
        }
        if(atom->position.z() < 0){
            atom->position[2] += m_systemSize[2];
            atom->initial_position[2] += m_systemSize[2];
        }
    }
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.

    vec3 totalMomentum(0,0,0);

    for (Atom* atom : m_atoms) {
        totalMomentum += atom->velocity * atom->mass();
    }
    totalMomentum /= m_atoms.size();

    for (Atom* atom : m_atoms) {
        atom->velocity -= totalMomentum / atom->mass();
    }

}

void System::resetForce(){
    for(Atom *atom : m_atoms){
        atom->resetForce();
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    //Random::randomSeed();

    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).
    for(int i = 0; i < numberOfUnitCellsEachDimension; i++)
        for(int j = 0; j < numberOfUnitCellsEachDimension; j++)
            for(int k = 0; k < numberOfUnitCellsEachDimension; k++)
                for(int l = 0; l < 4;  l++)
                {
                    double x,y,z;
                    Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                    if(l == 0){
                        x = i * latticeConstant;
                        y = j * latticeConstant;
                        z = k * latticeConstant;
                    }
                    else if(l == 1){
                        x = i * latticeConstant+latticeConstant/2.;
                        y = j * latticeConstant+latticeConstant/2.;
                        z = k * latticeConstant;
                    }
                    else if(l == 2){
                        x = i * latticeConstant+latticeConstant/2.;
                        y = j * latticeConstant;
                        z = k * latticeConstant+latticeConstant/2.;
                    }
                    else if(l == 3){
                        x = i * latticeConstant;
                        y = j * latticeConstant+latticeConstant/2.;
                        z = k * latticeConstant+latticeConstant/2.;
                    }

                    atom->position.set(x,y,z);
                    atom->initial_position.set(x,y,z);

                    double boltzmannConstant = 1.0; // In these units, the boltzmann constant equals 1
                    double standardDeviation = sqrt(boltzmannConstant*temperature/atom->mass());
                    double vx = Random::nextGaussian(0, standardDeviation);
                    double vy = Random::nextGaussian(0, standardDeviation);
                    double vz = Random::nextGaussian(0, standardDeviation);
                    atom->velocity.set(vx, vy, vz);
                    m_atoms.push_back(atom);
                }

    setSystemSize(vec3(latticeConstant*numberOfUnitCellsEachDimension,
                       latticeConstant*numberOfUnitCellsEachDimension,
                       latticeConstant*numberOfUnitCellsEachDimension));
}




void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();

    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}


void System::step(double dt) {
    m_integrator.integrate(*this, dt);
    m_steps++;
    m_time += dt;
}

