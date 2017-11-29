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
 * \param positions Positions of the particles
 * \param n1 Number of cells in the first direction
 * \param n2 Number of cells in the second direction
 */
void visuThread(std::shared_ptr<const PositionVec> positions,
                const long n1, const long n2) {
	// Initializations
    sf::RenderWindow window;
	float scale;
	int windowWidth, windowHeight;
	calcScale(scale, windowWidth, windowHeight, n1, n2);
    window.create(sf::VideoMode(windowWidth, windowHeight),
	              "Force in crystal");

    sf::CircleShape circle(Visu::circleRad);
	float a = scale * n1 * Hex::ux;
	sf::ConvexShape boundary(4);
	boundary.setPoint(0, sf::Vector2f(0, 0));
	boundary.setPoint(1, sf::Vector2f(a, 0));
	boundary.setPoint(2, sf::Vector2f(windowWidth, windowHeight));
	boundary.setPoint(3, sf::Vector2f(windowWidth - a, windowHeight));
	boundary.setFillColor(sf::Color::Yellow);

    window.setFramerateLimit(Visu::FPS);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

		// Update
        window.clear(sf::Color::White);
		window.draw(boundary);

		// Particle 0
		double x = (*positions)[0][0] * scale + windowWidth / 2; 
		double y = (*positions)[0][1] * scale + windowHeight / 2; 
		pbcHex(x, y, n1 * scale, n2 * scale);
		circle.setFillColor(sf::Color::Red);
		circle.setPosition(x, y);
		window.draw(circle);
		circle.setFillColor(sf::Color::Blue);

		for (long i = 1 ; i < n1 * n2 ; ++i) {
			double x = (*positions)[i][0] * scale + windowWidth / 2; 
			double y = (*positions)[i][1] * scale + windowHeight / 2; 
			pbcHex(x, y, n1 * scale, n2 * scale);
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
 * \param scale Scale for physical space to pixels
 * \param windowWidth Width of the window
 * \param windowHeight Height of the window
 * \param n1 Number of cells in the first direction
 * \param n2 Number of cells in the second direction
 */
void calcScale(float &scale, int &windowWidth, int &windowHeight,
               const long n1, const long n2) {
	float scale1 = Visu::windowSizeMax / (n1 * Hex::ux +  n2 * Hex::vx);
	float scale2 = Visu::windowSizeMax / (n1 * Hex::uy +  n2 * Hex::vy);
	scale = (scale1 < scale2) ? scale1 : scale2; // Minimum
	windowWidth = scale * (n1 * Hex::ux +  n2 * Hex::vx);
	windowHeight = scale * (n1 * Hex::uy +  n2 * Hex::vy);
}
