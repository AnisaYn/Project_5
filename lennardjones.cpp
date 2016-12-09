#include "lennardjones.h"
#include "system.h"
#include <cmath>

double LennardJones::potentialEnergy() const
{
    return m_potentialEnergy;
}

double LennardJones::sigma() const
{
    return m_sigma;
}

void LennardJones::setSigma(double sigma)
{
    m_sigma = sigma;
}

double LennardJones::epsilon() const
{
    return m_epsilon;
}

void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void LennardJones::calculateForces(System &system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop
    double sigma6 = pow(m_sigma, 6);
    double sigma12 = pow(m_sigma, 12);
    double L = system.systemSize().x();

    for(int i = 0; i < system.atoms().size(); i++)
    {
        Atom *atomi = system.atoms()[i];
        for(int j = i+1; j< system.atoms().size(); j++)
        {
            Atom *atomj = system.atoms()[j];
            double dx = atomi->position[0] - atomj->position[0];
            double dy = atomi->position[1] - atomj->position[1];
            double dz = atomi->position[2] - atomj->position[2];

            // Applying minimum image convention
            if(dx > 0.5*L) dx -= L;
            else if(dx < -0.5*L) dx += L;
            if(dy > 0.5*L) dy -= L;
            else if(dy < -0.5*L) dy += L;
            if(dz > 0.5*L) dz -= L;
            else if(dz < -0.5*L) dz += L;

            double r2 = dx*dx + dy*dy + dz*dz;
            double r2Inv = 1.0 / r2;
            double r6Inv = r2Inv*r2Inv*r2Inv;
            double r12Inv = r6Inv*r6Inv;
            double force = -24*m_epsilon*(sigma6*r6Inv - 2*sigma12*r12Inv)*r2Inv;
            atomi->force[0] += force*dx;
            atomi->force[1] += force*dy;
            atomi->force[2] += force*dz;

            atomj->force[0] -= force*dx;
            atomj->force[1] -= force*dy;
            atomj->force[2] -= force*dz;

            m_potentialEnergy += 4*m_epsilon*(sigma12*r12Inv - sigma6*r6Inv);
        }
    }
}
