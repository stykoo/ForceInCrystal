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
 * \brief Simulation of the system
 *
 * Implementation of the methods of the class Simul to simulate
 * a 2d hexagonal lattice of Brownian particles with one of them
 * submitted to an external force.
*/

#include <exception>
#include <thread>
#include <boost/program_options.hpp>
#include "simul.h"
#include "state.h"
#include "visu.h"

namespace po = boost::program_options;

/*!
 * \brief Constructor of Simul
 *
 * Initializes the parameters of structure Simul
 * from the command-line arguments using boost::program_options.
 *
 * \param argc Number of arguments
 * \param argv Arguments
 */
Simul::Simul(int argc, char **argv) {
	status = SIMUL_INIT_SUCCESS;

	po::options_description opts("Options");
	opts.add_options()
		("n1", po::value<long>(&n1)->required(),
		 "Number of cells in the first direction")
		("n2", po::value<long>(&n2)->required(),
		 "Number of cells in the first direction")
		("eps", po::value<double>(&potStrength)->required(),
		 "Strength of the potential")
		("T", po::value<double>(&temperature)->required(), "Temperature")
		("N", po::value<long>(&nbIters)->required(),
		 "Number of time iterations")
		("dt", po::value<double>(&dt)->required(), "Timestep")
		("sleep", po::value<int>(&sleep)->default_value(0),
		 "Number of milliseconds to sleep for between iterations")
		("help,h", "Print help message and exit")
		;

	try {
		po::variables_map vars;
		po::store(po::parse_command_line(argc, argv, opts), vars);

		// Display help and exit
		if (vars.count("help")) {
			std::cout << "Usage: " << argv[0] << " options\n";
			std::cout << opts << std::endl;
			status = SIMUL_INIT_HELP;
			return;
		}

        po::notify(vars);
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		status = SIMUL_INIT_FAILED;
		return;
	}

	// Check if the values of the parameters are allowed
	if (notPositive(n1, "n1") || notPositive(n2, "n2")
		|| notPositive(potStrength, "eps") || notPositive(temperature, "T")
		|| notPositive(dt, "dt") || notPositive(nbIters, "N")) {
		status = SIMUL_INIT_FAILED;
		return;
	}
}

//! Run the simulation
void Simul::run() {
	if (status != SIMUL_INIT_SUCCESS) {
		std::cerr << "You should not be runing a failed simulation..."
		          << std::endl;
		return;
	}

	// Initialize the state of the system
	State state(n1, n2, potStrength, temperature, dt);
	std::shared_ptr<const PositionVec> positions = state.getPositions();
	
	// Start thread for visualization
	std::thread thVisu(visuThread, positions, n1, n2); 

	// Time evolution
	for (long t = 0 ; t < nbIters ; ++t) {
		state.evolve();
		if (sleep > 0) {
			std::this_thread::sleep_for(std::chrono::milliseconds(sleep));
		}
	}

	thVisu.join();
}
