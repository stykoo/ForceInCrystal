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
	         const double _temperature, const double _dt) :
	n1(_n1), n2(_n2), potStrength(_potStrength), dt(_dt),
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
			(*positions)[ind][0] = i * Hex::ux + j * Hex::vx;
			(*positions)[ind][1] = i * Hex::uy + j * Hex::vy;
		}
	}
}

/*!
 * \brief Do one time step
 *
 * Evolve the system for one time step according to coupled Langevin equation.
 */
void State::evolve() {
	// TODO
}
