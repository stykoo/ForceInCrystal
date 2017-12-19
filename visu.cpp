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
 * \brief Visualization of the system
 *
 * The system is visualized using the SFML library.
*/

#include <SFML/Graphics.hpp>
#include "visu.h"

/*!
 * \brief Thread for visualization.
 *
 * Open a window, draw the particles and update their
 * positions at a certain number of FPS while the simulation is runing.
 *
 * \param state State of the system
 * \param n1 Number of cells in the first direction
 * \param n2 Number of cells in the second direction
 */
void visuThread(const State *state, const long n1, const long n2) {
	double Lx = n1;
	double Ly = n2 * Hex::vy;

	// Initializations
    sf::RenderWindow window;
	float scale = calcScale(Lx, Ly);
	int windowWidth = Lx * scale;
	int windowHeight = Ly * scale;
    window.create(sf::VideoMode(windowWidth, windowHeight),
	              "Force in crystal");

    sf::CircleShape circle(Visu::circleRad);

    window.setFramerateLimit(Visu::FPS);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

		// Update
        window.clear(sf::Color::Yellow);

		// Particle 0
		double x = state->getPos(0)[0] * scale + windowWidth / 2; 
		pbc(x, windowWidth);
		double y = state->getPos(0)[1] * scale + windowHeight / 2; 
		pbc(y, windowHeight);
		circle.setFillColor(sf::Color::Red);
		circle.setPosition(x, y);
		window.draw(circle);
		circle.setFillColor(sf::Color::Blue);

		for (long i = 1 ; i < n1 * n2 ; ++i) {
			double x = state->getPos(i)[0] * scale + windowWidth / 2; 
			pbc(x, windowWidth);
			double y = state->getPos(i)[1] * scale + windowHeight / 2; 
			pbc(y, windowHeight);
			circle.setPosition(x, y);
            window.draw(circle);
		}
        window.display();
    }
}

/*!
 * \brief Compute the scale to go from physical unit to screen pixel unit.
 *
 * Compute the scale to go from physical unit to screen pixel unit
 * and the width and height of the window.
 *
 * \param Lx Length in the x direction
 * \param Ly Length in the y direction
 */
float calcScale(const double Lx, const double Ly) {
	float scale1 = Visu::windowSizeMax / Lx;
	float scale2 = Visu::windowSizeMax / Ly;
	return (scale1 < scale2) ? scale1 : scale2; // Minimum
}

