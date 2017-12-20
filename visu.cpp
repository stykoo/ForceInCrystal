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

#include "visu.h"
#include <iostream>

Visu::Visu(const State *state_, const long n1_, const long n2_) :
	state(state_), n1(n1_), n2(n2_), Lx(n1_), Ly(n2_ * Hex::vy),
	scale(std::min(windowSizeMax / Lx, windowSizeMax / Ly)),
	windowWidth(Lx * scale), windowHeight(Ly * scale),
	shiftX(windowWidth / 2), shiftY(windowHeight / 2) {
}

/*!
 * \brief Thread for visualization.
 *
 * Open a window, draw the particles and update their
 * positions at a certain number of FPS while the simulation is runing.
 */
void Visu::run() {
    sf::RenderWindow window;
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

		// Voronoi tesselation
		sf::VertexArray lines(sf::Lines);
		genLinesVoronoi(lines);

		// Update
        window.clear(sf::Color::Yellow);

		// Lines of Voronoi tesselation
		window.draw(lines);

		// Particle 0
		double x = state->getPos(0)[0] * scale + shiftX; 
		pbc(x, windowWidth);
		double y = state->getPos(0)[1] * scale + shiftY; 
		pbc(y, windowHeight);
		circle.setFillColor(sf::Color::Red);
		circle.setPosition(x, y);
		window.draw(circle);
		circle.setFillColor(sf::Color::Blue);

		for (long i = 1 ; i < n1 * n2 ; ++i) {
			double x = state->getPos(i)[0] * scale + shiftX; 
			pbc(x, windowWidth);
			double y = state->getPos(i)[1] * scale + shiftY; 
			pbc(y, windowHeight);
			circle.setPosition(x, y);
            window.draw(circle);
		}
        window.display();
    }
}

void Visu::genLinesVoronoi(sf::VertexArray &lines) {
	VoronoiList vor = state->computeVoronoi();
	double x0 = state->getPos(0)[0];
	double y0 = state->getPos(0)[1];

	for (CGAL_K::Segment_2 s : vor) {
		double x1 = scale * (s[0][0] + x0) + shiftX;
		double y1 = scale * (s[0][1] + y0) + shiftY;
		double x2 = scale * (s[1][0] + x0) + shiftX;
		double y2 = scale * (s[1][1] + y0) + shiftY;
		int ox1 = 0, oy1 = 0, ox2 = 0, oy2 = 0;
		pbcOffset(x1, ox1, windowWidth);
		pbcOffset(y1, oy1, windowHeight);
		pbcOffset(x2, ox2, windowWidth);
		pbcOffset(y2, oy2, windowHeight);

		if (ox1 == ox2 && oy1 == oy2) {
				lines.append(sf::Vertex(sf::Vector2f(x1, y1),
				                        sf::Color::Black));
				lines.append(sf::Vertex(sf::Vector2f(x2, y2),
				                        sf::Color::Black));
		} 
		// We choose not to plot the segments crossing the boundaries
		// of the window. We could also plot two segments for these.
	}
}
