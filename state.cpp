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
 * \param _temperature Temperature
 * \param _force External force on particle 0
 * \param _dt Timestep
 * \param _screening Screening length
 */
State::State(const long _n1, const long _n2,
	         const double _temperature, const double _force, const double _dt,
			 const double _screening) :
	n1(_n1), n2(_n2), fx(_force), fy(0), dt(_dt), screening(_screening),
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
		// External force on particle 0
		(*positions)[0][0] += dt * fx;
		(*positions)[0][1] += dt * fy;
		pbcHex((*positions)[i][0], (*positions)[i][1], n1, n2);
	}
}

/* \brief Compute the forces between the particles.
 *
 * Implement a screened dipole-dipole interaction.
 */
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
			double u = (3.0 + dr / screening)
				        * std::exp(- dr / screening) / (dr2 * dr2);
			double fx = u * dx / dr;
			double fy = u * dy / dr;

			forces[i][0] += fx;
			forces[j][0] -= fx;
			forces[i][1] += fy;
			forces[j][1] -= fy;
        }
    }
}

/*! 
 * \brief Periodic boundary conditions on a segment
 * 
 * Update x to be between 0 and L.
 *
 * \param x Value
 * \param L Length of the box
 */
void pbc(double &x, const double L) {
	x -= L * std::floor(x / L);
}

/*! 
 * \brief Periodic boundary conditions on a segment (symmetric)
 * 
 * Update x to be between -L/2 and L/2.
 *
 * \param x Value
 * \param L Length of the box
 */
void pbcSym(double &x, const double L) {
	x -= L * std::round(x / L);
}

/*! 
 * \brief Periodic boundary conditions on a hexagonal lattice
 * 
 * Update x to be between x and y to be in the cell of lengths L1 and L2.
 *
 * \param x Value 1
 * \param y Value 2
 * \param L1 Length of the cell in the first direction
 * \param L2 Length of the cell in the first direction
 */
void pbcHex(double &x, double &y, const double L1, const double L2) {
	double a = x * Hex::inv11 + y * Hex::inv12;	
	double b = x * Hex::inv21 + y * Hex::inv22;	
	pbc(a, L1);
	pbc(b, L2);
	x = a * Hex::ux + b * Hex::vx;
	y = a * Hex::uy + b * Hex::vy;
}

/*! 
 * \brief Periodic boundary conditions on a hexagonal lattice (symmetric)
 * 
 * Update x to be between x and y to be in the cell of lengths L1 and L2
 * with (0, 0) being the center of the cell.
 *
 * \param x Value 1
 * \param y Value 2
 * \param L1 Length of the cell in the first direction
 * \param L2 Length of the cell in the first direction
 */
void pbcHexSym(double &x, double &y, const double L1, const double L2) {
	double a = x * Hex::inv11 + y * Hex::inv12;	
	double b = x * Hex::inv21 + y * Hex::inv22;	
	pbcSym(a, L1);
	pbcSym(b, L2);
	x = a * Hex::ux + b * Hex::vx;
	y = a * Hex::uy + b * Hex::vy;
}
