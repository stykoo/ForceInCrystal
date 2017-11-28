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
 * \file simul.h
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \version 1.0
 * \date 2017-11-27
 * \brief Simulation of the system
 *
 * Header file for simul.cpp.
 * It defines the class Simul.
 */

#ifndef FORCEINCRYSTAL_SIMUL_H_
#define FORCEINCRYSTAL_SIMUL_H_

#include <string>
#include <iostream>

//! State of the simulation after initialization
enum SimulInitStatus {
	SIMUL_INIT_SUCCESS, //!< Successful initialization
	SIMUL_INIT_HELP, //!< Display help
	SIMUL_INIT_FAILED //!< Failed initialization
};

/*!
 * \brief Class for simulation
 *
 * This class takes care of both the initialization
 * and the time evolution of the system.
 */
class Simul {
	public:
		Simul(int argc, char **argv); //!< Constructor from arguments
		void run(); //!< Run the simulation

		//! Get initialization status
		SimulInitStatus getStatus() const { return status; }

	private:
		long n1; //!< Number of cells in the first direction
		long n2; //!< Number of cells in the second direction
		double potStrength; //!< Strength of the potential
		double temperature; //!< Temperature
		double force; //!< External force on particle 0
		double dt; //!< Timestep
		long nbIters; //! Number of time iterations
		double screening; //! Screening length

		int sleep; //!< Number of milliseconds to sleep for between iterations
		SimulInitStatus status; //!< Status after initialization
};

/*!
 * \brief Check if variable is positive.
 *
 * Returns true and prints an error message if the variable
 * is not strictly positive.
 *
 * \param a Variable
 * \param std::string Name of variable
 * \return false if the variable is strictly positive, true otherwise
 */
template<typename T>
bool notPositive(const T &a, std::string name) {
	if (a <= T(0)) {
		std::cerr << "Error: " << name << " should be strictly positive."
		          << std::endl;
		return true;
	}
	return false;
}

#endif // FORCEINCRYSTAL_SIMUL_H_
