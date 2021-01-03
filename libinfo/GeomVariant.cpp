// Vagabond
// Copyright (C) 2019 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#include "GeomVariantDefs.h"
#include <cmath>
#include <iostream>

double GeomVariant::getAngle(double *table, double phi, double psi)
{
	/* get these -180 to 170 values to vary between 0 and 35 */
	phi = (phi + 180) / 10;
	psi = (psi + 180) / 10;

	/* psi is the fast axis */
	int lookup = (int)phi * 36 + (int)psi;

	return table[lookup];
}

double GeomVariant::interpolateAngle(double *table, double phi, double psi)
{
	/* get these -180 to 170 values to vary between 0 and 35 */
	phi = (phi + 180) / 10;
	psi = (psi + 180) / 10;
	
	double phi1 = phi + 1;
	double psi1 = psi + 1;
	
	if (phi >= 36) phi -= 36;
	if (psi >= 36) psi -= 36;
	if (phi1 >= 36) phi1 -= 36;
	if (psi1 >= 36) psi1 -= 36;

	/* psi is the fast axis */
	int lookup00 = (int)phi * 36 + (int)psi;
	int lookup01 = (int)phi * 36 + (int)psi1;
	int lookup10 = (int)phi1 * 36 + (int)psi;
	int lookup11 = (int)phi1 * 36 + (int)psi1;
	
	double rphi = phi - lrint(phi);
	double rpsi = psi - lrint(psi);

	double result = table[lookup00];
	result += rphi * (table[lookup10] - table[lookup00]);
	result += rpsi * (table[lookup01] - table[lookup00]);
	result += rpsi * rphi * (table[lookup11] + table[lookup00]
	                         - table[lookup10] - table[lookup01]);
	
	return result;
}

