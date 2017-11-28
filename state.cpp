/*
Copyright (C) Université Pierre et Marie Curie (2017)
Contributor: Alexis Poncet <aponcet@lptmc.jussieu.fr>

This file is part of ForceInCrystal.

ForceInCrystal is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ForceInCrystal is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ForceInCrystal.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file simul.cpp
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \version 1.0
 * \date 2017-11-27
 * \brief State of the system
 *
 * Implementation of the methods of the class State to simulate
 * of a 2d hexagonal lattice of Brownian particles with one of them
 * submitted to an external force.
*/

#include <cmath>
#include <chrono>
#include <iostream>
#include "state.h"

/*!
 * \brief Constructor of State
 *
 * Initializes the state of the system
 * as particles in a 
 *
 * \param _n1 Number of cells in the first direction
 * \param _n2 Number of cells in the second direction
 * \param _potStrength Strength of the potential
 * \param _temperature Temperature
 * \param _dt Timestep
 */
State::State(const long _n1, const long _n2, const double _potStrength,
	         const double _temperature, const double _dt,
			 const double _screening) :
	n1(_n1), n2(_n2), potStrength(_potStrength), dt(_dt),
	screening(_screening),
	// We initialize the gaussian noise from the temperature
	gaussianNoise(0.0, std::sqrt(2.0 * _temperature * dt)),
	// We seed the RNG with the current time
	rng(std::chrono::system_clock::now().time_since_epoch().count())
{
	// Put the particles on a hexagonal lattice	
	positions.reset(new PositionVec(n1 * n2));
	for (long i = 0 ; i < n1 ; ++i) {
		for (long j = 0 ; j < n2 ; ++j) {
			long ind = i * n2 + j;
			(*positions)[ind][0] = (i + 0.5) * Hex::ux + (j + 0.5) * Hex::vx;
			(*positions)[ind][1] = (i + 0.5) * Hex::uy + (j + 0.5) * Hex::vy;
		}
	}

	forces.resize(n1 * n2);
}

/*!
 * \brief Do one time step
 *
 * Evolve the system for one time step according to coupled Langevin equation.
 */
void State::evolve() {
	calcInternalForces();
	for (long i = 0 ; i < n1 * n2 ; ++i) {
		// Internal forces + Gaussian noise
		(*positions)[i][0] += dt * forces[i][0] + gaussianNoise(rng);
		(*positions)[i][1] += dt * forces[i][1] + gaussianNoise(rng);
		pbcHex((*positions)[i][0], (*positions)[i][1], n1, n2);
	}
}

void State::calcInternalForces() {
    for (long i = 0 ; i < n1 * n2 ; ++i) {
		forces[i][0] = 0;
		forces[i][1] = 0;
    }

    for (long i = 0 ; i < n1 * n2 ; ++i) {
        for (long j = i + 1 ; j < n1 * n2 ; ++j) {
			double dx = (*positions)[i][0] - (*positions)[j][0];
			double dy = (*positions)[i][1] - (*positions)[j][1];
			// We want the value of dx to be between -n1/2 and n1/2
			pbcHexSym(dx, dy, n1, n2);
			double dr2 = dx * dx + dy * dy;
			double dr = std::sqrt(dr2);
			double t = - potStrength * (3.0 + dr / screening)
			           * std::exp(- dr / screening) / (dr2 * dr2 * dr);
			forces[i][0] += t * dx;
			forces[j][0] -= t * dx;
			forces[i][1] += t * dy;
			forces[j][1] -= t * dy;
        }
    }
	/*
	for (long i = 0 ; i < n1 ; ++i) {
		for (long j = 0 ; j < n2 ; ++j) {
			long ind = i * n2 + j;
			std::cout << forces[ind][0]  << " " << forces[ind][1] << "\n";
		}
		std::cout << "\n";
    }
	*/
}

void pbc(double &x, const double L) {
	x -= L * std::floor(x / L);
}

void pbcSym(double &x, const double L) {
	x -= L * std::round(x / L);
}

void pbcHex(double &x, double &y, const double L1, const double L2) {
	double a = x * Hex::inv11 + y * Hex::inv12;	
	double b = x * Hex::inv21 + y * Hex::inv22;	
	pbc(a, L1);
	pbc(b, L2);
	x = a * Hex::ux + b * Hex::vx;
	y = a * Hex::uy + b * Hex::vy;
}

void pbcHexSym(double &x, double &y, const double L1, const double L2) {
	double a = x * Hex::inv11 + y * Hex::inv12;	
	double b = x * Hex::inv21 + y * Hex::inv22;	
	pbcSym(a, L1);
	pbcSym(b, L2);
	x = a * Hex::ux + b * Hex::vx;
	y = a * Hex::uy + b * Hex::vy;
}
