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
 * \file state.h
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \version 1.0
 * \date 2017-11-27
 * \brief State of the system
 *
 * Header file for state.h.
 * It defines the class State.
 */

#ifndef FORCEINCRYSTAL_STATE_H_
#define FORCEINCRYSTAL_STATE_H_

#include <vector>
#include <array>
#include <random>
#include <memory>

//! Basis vectors of hexagonal lattice
namespace Hex {
	const double ux = 1.0; //!< 1st component of 1st vector
	const double uy = 0; //!< 2nd component of 1st vector
	const double vx = 0.5; //!< 1st component of 2nd vector
	const double vy = 0.86602540378443864676; //!< 2nd component of 2nd vector 
	const double inv11 = 1.0; //!< Element (1, 1) of the inverse matrix
	//! Element (1, 2) of the inverse matrix
	const double inv12 = -0.57735026918962576451;
	const double inv21 = 0; //!< Element (2, 1) of the inverted matrix
	//! Element (2, 2) of the inverse matrix
	const double inv22 =  1.15470053837925152902;
}

//! Type name for vector of positions
typedef std::vector< std::array<double, 2> > PositionVec;

/*!
 * \brief Class for the state of the system
 *
 * This class takes care of the initialization
 * and the evolution of the state of the system.
 */
class State {
	public:
		//! Constructor of State
		State(const long _n1, const long _n2, const double _potStrength,
		      const double _temperature, const double _force, const double _dt,
			  const double _screening);
		void evolve(); //!< Do one time step

		//! Return a pointer on the positions of the particles
		std::shared_ptr<const PositionVec> getPositions() const {
			return positions; 
		}

	private:
		void calcInternalForces(); //!< Compute internal forces

		const long n1; //!< Number of cells in the first direction
		const long n2; //!< Number of cells in the second direction
		const double potStrength; //!< Strength of the potential
		double fx, fy; //!< External force on particle 0
		const double dt; //!< Timestep
		const double screening; //!< Screening length

		std::normal_distribution<double> gaussianNoise;  //!< Gaussian noise
		std::mt19937 rng; //! Random number generator

		//! Positions of the particles
		// Shared pointer because it will also be used for visualization
		std::shared_ptr<PositionVec> positions;
		PositionVec forces; //!< Internal forces
};

//! Periodic boundary conditions on a segment
void pbc(double &x, const double L);
//! Periodic boundary conditions on a segment (symmetric)
void pbcSym(double &x, const double L);
//! Periodic boundary conditions on a hexagonal lattice
void pbcHex(double &x, double &y, const double L1, const double L2);
//! Periodic boundary conditions on a hexagonal lattice (symmetric)
void pbcHexSym(double &x, double &y, const double L1, const double L2);

#endif // FORCEINCRYSTAL_STATE_H_
