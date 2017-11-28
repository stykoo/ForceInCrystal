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
 * \date 2017-11-28
 * \brief Unit testing for state.cpp
*/

#define BOOST_TEST_MODULE My Test
#include <boost/test/included/unit_test.hpp>
#include "state.h"

//! Test of pbc(x, L)
BOOST_AUTO_TEST_CASE(pbc_test)
{
	double a = 0.105, x = 2.105, y = 10.105, z = -1.895;
	pbc(a, 2.0);
	pbc(x, 2.0);
	pbc(y, 2.0);
	pbc(z, 2.0);
	BOOST_CHECK_CLOSE(a, 0.105, 1e-5);
	BOOST_CHECK_CLOSE(x, 0.105, 1e-5);
	BOOST_CHECK_CLOSE(y, 0.105, 1e-5);
	BOOST_CHECK_CLOSE(z, 0.105, 1e-5);
}

//! Test of pbcHex(x, y, L1, L2)
BOOST_AUTO_TEST_CASE(pbc_hex_test)
{
	double L1 = 3.0, L2 = 2.0;
	double x, y;

	x = 2.0; y = 1.0;
	pbcHex(x, y, L1, L2);
	BOOST_CHECK_CLOSE(x, 2.0, 1e-5); BOOST_CHECK_CLOSE(y, 1.0, 1e-5);

	x = 4.0; y = 1.0;
	pbcHex(x, y, L1, L2);
	BOOST_CHECK_CLOSE(x, 1.0, 1e-5); BOOST_CHECK_CLOSE(y, 1.0, 1e-5);

	x = -5.0; y = 1.0;
	pbcHex(x, y, L1, L2);
	BOOST_CHECK_CLOSE(x, 1.0, 1e-5); BOOST_CHECK_CLOSE(y, 1.0, 1e-5);

	x = 2.0; y = 2.0;
	pbcHex(x, y, L1, L2);
	BOOST_CHECK_CLOSE(x, 2.0 - L2 * 0.5, 1e-5);
	BOOST_CHECK_CLOSE(y, 2.0 - L2 * 0.86602540378443864676, 1e-5);

	x = 1.0; y = -1.0;
	pbcHex(x, y, L1, L2);
	BOOST_CHECK_CLOSE(x, 1.0 + L2 * 0.5, 1e-5);
	BOOST_CHECK_CLOSE(y, -1.0 + L2 * 0.86602540378443864676, 1e-5);
}
