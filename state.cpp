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
 * \param _temperature Temperature
 * \param _fv External force or velocity on particle 0
 * \param _angle Angle in degrees for external force or velocity on particle 0
 * \param _dt Timestep
 * \param _screening Screening length
 * \param _evolType Constant force or velocity
 */
State::State(const long _n1, const long _n2,
	         const double _temperature, const double _fv, const double _angle,
			 const double _dt, const double _screening,
			 const StateEvolType _evolType, const StatePBCType _pbcType) :
	n1(_n1), n2(_n2),
	// Caution: the length in y in not n2!
	Lx(_n1), Ly(_n2 * Hex::vy),
	fvx(_fv * std::cos(_angle * M_PI / 180)),
	fvy(_fv * std::sin(_angle * M_PI / 180)),
	dt(_dt), screening(_screening), evolType(_evolType), pbcType(_pbcType),
	// We initialize the gaussian noise from the temperature
	gaussianNoise(0.0, std::sqrt(2.0 * _temperature * dt)),
	// We seed the RNG with the current time
	rng(std::chrono::system_clock::now().time_since_epoch().count())
{
	positions = std::make_shared<PositionVec>();
	(*positions)[0].resize(n1 * n2);
	(*positions)[1].resize(n1 * n2);

	// Put the particles on a hexagonal lattice	
	if (pbcType == SQUARE_PBC) {
		for (long i = 0 ; i < n1 ; ++i) {
			for (long j = 0 ; j < n2 ; ++j) {
				long ind = i * n2 + j;
				(*positions)[0][ind] = i + 0.5 * (j % 2) + 0.25;
				(*positions)[1][ind] = (j + 0.5) * Hex::vy;
			}
		}
	} else if (pbcType == HEX_PBC) {
		for (long i = 0 ; i < n1 ; ++i) {
			for (long j = 0 ; j < n2 ; ++j) {
				long ind = i * n2 + j;
				(*positions)[0][ind] = (i + 0.5) * Hex::ux
				                       + (j + 0.5) * Hex::vx;
				(*positions)[1][ind] = (i + 0.5) * Hex::uy
				                       + (j + 0.5) * Hex::vy;
			}
		}
	}
	enforcePBC();

	forces[0].resize(n1 * n2);
	forces[1].resize(n1 * n2);
}

/*!
 * \brief Do one time step
 *
 * Evolve the system for one time step according to coupled Langevin equation.
 */
void State::evolve() {
	calcInternalForces();
	// Particles other than 0
	for (long i = 1 ; i < n1 * n2 ; ++i) {
		// Internal forces + Gaussian noise
		(*positions)[0][i] += dt * forces[0][i] + gaussianNoise(rng);
		(*positions)[1][i] += dt * forces[1][i] + gaussianNoise(rng);
	}

	// Particle 0
	if (evolType == CONSTANT_FORCE) {
		(*positions)[0][0] += dt * (forces[0][0] + fvx) + gaussianNoise(rng);
		(*positions)[1][0] += dt * (forces[1][0] + fvy) + gaussianNoise(rng);
	} else {
		(*positions)[0][0] += dt * fvx;
		(*positions)[1][0] += dt * fvy;
	}

	enforcePBC();
}

/* \brief Compute the forces between the particles.
 *
 * Implement a screened dipole-dipole interaction.
 */
void State::calcInternalForces() {
    for (long i = 0 ; i < n1 * n2 ; ++i) {
		forces[0][i] = 0;
		forces[1][i] = 0;
    }

    for (long i = 0 ; i < n1 * n2 ; ++i) {
        for (long j = i + 1 ; j < n1 * n2 ; ++j) {
			double dx = (*positions)[0][i] - (*positions)[0][j];
			double dy = (*positions)[1][i] - (*positions)[1][j];
			// We want the periodized interval to be centered in 0
			if (pbcType == SQUARE_PBC) {
				pbcSym(dx, Lx);
				pbcSym(dy, Ly);
			} else if (pbcType == HEX_PBC) {
				pbcHexSym(dx, dy, n1, n2);
			}
			double dr2 = dx * dx + dy * dy;
			double dr = std::sqrt(dr2);
			double u = (3.0 + dr / screening)
				        * std::exp(- dr / screening) / (dr2 * dr2);
			double fx = u * dx / dr;
			double fy = u * dy / dr;

			forces[0][i] += fx;
			forces[0][j] -= fx;
			forces[1][i] += fy;
			forces[1][j] -= fy;
        }
    }
}

/* \brief Enforce periodic boundary conditions
 *
 * Put the positions in an square or hexagonal box.
 */
void State::enforcePBC() {
	if (pbcType == SQUARE_PBC) {
		for (long i = 0 ; i < n1 * n2 ; ++i) {
			pbc((*positions)[0][i], Lx);
			pbc((*positions)[1][i], Ly);
		}
	} else if (pbcType == HEX_PBC) {
		for (long i = 0 ; i < n1 * n2 ; ++i) {
			pbcHex((*positions)[0][i], (*positions)[1][i], n1, n2);
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
