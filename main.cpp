/*
Copyright (C) Université Pierre et Marie Curie (2017)
Contributor: Alexis Poncet <aponcet@lptmc.jussieu.fr>

This software (ForceInCrystal) intends to investigate the behavior
of a 2d hexagonal lattice of Brownian particles when one of them
is submitted to an external force.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \mainpage Force in crystal
 *
 * ForceInCrystal is a software to investigate the behavior of a 2d
 * hexagonal lattice of Brownian particles when one of them is submitted
 * to an external force.
 */
/*!
 * \file main.cpp
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \version 1.0
 * \date 2017-11-27
 * \brief Main file of ForceInCrystal
*/

#include "simul.h"

/*!
 * \brief Main function
 *
 * Creates and run the simulation.
 */
int main(int argc, char **argv) {
	Simul simulation(argc, argv);

	if (simulation.getStatus() == SIMUL_INIT_HELP) {
		return 0;
	} else if (simulation.getStatus() == SIMUL_INIT_FAILED) {
		return 1;
	}

	simulation.print();
	simulation.run();

	return 0;
}
