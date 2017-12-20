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
			 const StateEvolType _evolType) :
	n1(_n1), n2(_n2),
	// Caution: the length in y in not n2!
	Lx(_n1), Ly(_n2 * Hex::vy),
	fvx(_fv * std::cos(_angle * M_PI / 180)),
	fvy(_fv * std::sin(_angle * M_PI / 180)),
	dt(_dt), screening(_screening), evolType(_evolType),
	// We initialize the gaussian noise from the temperature
	gaussianNoise(0.0, std::sqrt(2.0 * _temperature * dt)),
	// We seed the RNG with the current time
	rng(std::chrono::system_clock::now().time_since_epoch().count())
{
	positions.resize(n1 * n2);

	for (long i = 0 ; i < n1 ; ++i) {
		for (long j = 0 ; j < n2 ; ++j) {
			long ind = i * n2 + j;
			positions[ind][0] = i + 0.5 * (j % 2) + 0.25;
			positions[ind][1] = (j + 0.5) * Hex::vy;
		}
	}
	enforcePBC();

	forces.resize(n1 * n2);
}

/*!
 * \brief Return the Delaunay triangulation of the state
 *
 * The box is not considered as periodic: we center the triangulation
 * on particle 0.
 */
CGAL_DT State::computeDelaunay() const {
	std::list<CGAL_DT::Point> L;
	for (long i = 0 ; i < n1 * n2 ; ++i) {
		double x = positions[i][0] - positions[0][0];
		double y = positions[i][1] - positions[0][1];
		pbcSym(x, Lx);
		pbcSym(y, Ly);
		L.push_front(CGAL_DT::Point(x, y));
	}

	return CGAL_DT(L.begin(), L.end()); // Delaunay triangulation
}

/*!
 * \brief Return the Voronoi tesselation of the state
 *
 * The box is not considered as periodic: we center the tesselation
 * on particle 0.
 */
VoronoiList State::computeVoronoi() const {
	CGAL_DT T = computeDelaunay();
	// CGAL_K::Iso_rectangle_2 bbox(-Lx / 2, -Ly / 2, Lx / 2, Ly / 2);
	VoronoiList voronoi;

    for (CGAL_DT::Edge_iterator eit = T.edges_begin() ; eit != T.edges_end() ;
	     ++eit) {
		CGAL::Object obj = T.dual(eit);
		/*CGAL_K::Object inter = CGAL::intersection(bbox, obj);
		const CGAL_K::Segment_2* s =
			CGAL::object_cast<CGAL_K::Segment_2>(&inter); */
		const CGAL_K::Segment_2* s =
			CGAL::object_cast<CGAL_K::Segment_2>(&obj);
		if (s) {
			voronoi.push_back(*s);
		}
	}

	return voronoi; // This copy should be avoided
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
		positions[i][0] += dt * forces[i][0] + gaussianNoise(rng);
		positions[i][1] += dt * forces[i][1] + gaussianNoise(rng);
	}

	// Particle 0
	if (evolType == CONSTANT_FORCE) {
		positions[0][0] += dt * (forces[0][0] + fvx) + gaussianNoise(rng);
		positions[0][1] += dt * (forces[0][1] + fvy) + gaussianNoise(rng);
	} else {
		positions[0][0] += dt * fvx;
		positions[0][1] += dt * fvy;
	}

	enforcePBC();
}

/* \brief Compute the forces between the particles.
 *
 * Implement a screened dipole-dipole interaction.
 */
void State::calcInternalForces() {
    for (long i = 0 ; i < n1 * n2 ; ++i) {
		forces[i][0] = 0;
		forces[i][1] = 0;
    }

    for (long i = 0 ; i < n1 * n2 ; ++i) {
        for (long j = i + 1 ; j < n1 * n2 ; ++j) {
			double dx = positions[i][0] - positions[j][0];
			double dy = positions[i][1] - positions[j][1];
			// We want the periodized interval to be centered in 0
			pbcSym(dx, Lx);
			pbcSym(dy, Ly);
			double dr2 = dx * dx + dy * dy;
			double dr = std::sqrt(dr2);
			double u = (3.0 + dr / screening)
				        * std::exp(- dr / screening) / (dr2 * dr2);
			double fx = u * dx / dr;
			double fy = u * dy / dr;

			forces[i][0] += fx;
			forces[j][0] -= fx;
			forces[i][1] += fy;
			forces[j][1] -= fy;
        }
    }
}

/* \brief Enforce periodic boundary conditions
 *
 * Put the positions in a square box.
 */
void State::enforcePBC() {
	for (long i = 0 ; i < n1 * n2 ; ++i) {
		pbc(positions[i][0], Lx);
		pbc(positions[i][1], Ly);
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
