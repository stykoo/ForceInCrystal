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
 * \file visu.h
 * \author Alexis Poncet <aponcet@lptmc.jussieu.fr>
 * \version 1.0
 * \date 2017-11-27
 * \brief Visualization of the system
 *
 * Header file for visu.cpp.
*/

#ifndef FORCEINCRYSTAL_VISU_H_
#define FORCEINCRYSTAL_VISU_H_

#include <SFML/Graphics.hpp>
#include "state.h"

class Visu {
	public:
		Visu(const State *state, const long n1, const long n2);
		void run();

	private:
		//! Lines from Voronoi tesselation
		void genLinesVoronoi(sf::VertexArray &lines);

		// Constant variables for visualization
		const int windowSizeMax = 800; //!< Maximum size of the window
		const float circleRad = 5.0; //!< Radius of the particles on the screen
		const int FPS = 24; //!< Number of frames per second

		const State *state; //!< Pointer to the state of the system
		const long n1; //!< Number of cells in the first direction
		const long n2; //!< Number of cells in the second direction
		const double Lx; //!< Length in x
		const double Ly; //!< Length in y
		const float scale; //!< Scale from data to window
		const int windowWidth; //< Width of the window
		const int windowHeight; //< Height of the window
		const int shiftX; //< Shift in the x direction
		const int shiftY; //< Shift in the y direction
};

#endif // FORCEINCRYSTAL_VISU_H_
