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
 * Header file for simul.cpp.
*/

#ifndef FORCEINCRYSTAL_VISU_H_
#define FORCEINCRYSTAL_VISU_H_

#include "state.h"

namespace Visu {
	const int windowSizeMax = 900;
	const float circleRad = 5.0;
	const int FPS = 24;
}

//! Thread for visualizing the particles
void visuThread(std::shared_ptr<const PositionVec> positions,
                const long n1, const long n2); 

//! Compute scale from number of particles
void calcScale(float &scale, int &windowWidth, int &windowHeight,
               const long n1, const long n2);

#endif // FORCEINCRYSTAL_VISU_H_
