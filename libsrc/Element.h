//
//  Element.h
//  vagabond
//
//  Created by Helen Ginn on 16/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Element__
#define __vagabond__Element__

#include <stdio.h>
#include <string>
#include <vector>
#include "shared_ptrs.h"
#include "Distributor.h"

/**
 * \class Element
 * \brief Maintains the scattering factors and pre-calculated distributions
 * for a single element of the periodic table.
 * 
 * As long as Element::setupElements() has been called, element objects
 * can be accessed using their capitalised two letter symbol as an argument
 * to Element::getElement().
 */

class Element : public Distributor
{
public:

	/** To be called once at the beginning: set up all known elements which
	 * can then be accessed using the static function getElement() */
	static void setupElements();

	~Element() {}
	
	/** Get the appropriate Element object for a given capitalised two-letter
	  * symbol as seen on the periodic table */
	static ElementPtr getElement(std::string symbol);
	
	/** Get the FFT distribution for this atom, using the scattering
	 * 	factors from International Tables F */
	virtual FFTPtr getDistribution(bool = false, int new_n = -1);
	FFTPtr getMask();

	/** Returns one- or two-letter abbreviation of element as on the periodic
	  * table, fully capitalised */
	std::string getSymbol()
	{
		return _symbol;
	}

	/** Lower-case fully written out element name, e.g. "nitrogen" */
	std::string getName()
	{
		return _name;
	}

	/** Number of electrons for an atom of a given element */
	double electronCount()
	{
		return _electrons;
	}

	/** Returns "Element" for the class name. */
	virtual std::string getClassName()
	{
		return "Element";
	}

	/** Pass in a std::vector of atoms and receive an std::vector of unique
	 * Elements covering all atoms in the group. */
	static std::vector<ElementPtr> elementList(std::vector<AtomPtr> atoms);
protected:
	/** Return an atomic scattering factor value for a given position
	 * @param obj pointer to Element object which has been cast to void
	 * @param x x position in inverse Angstroms
	 * @param y y position in inverse Angstroms
	 * @param z z position in inverse Angstroms
	 */
	static double getVoxelValue(void *obj, double x, double y, double z);
	static double getSolventMaskValue(void *obj, double x, double y, double z);

private:
	Element(std::string symbol, std::string name, double electrons,
	        const float *scatter);
	std::string _symbol;
	std::string _name;
	double _electrons;

	float _scattering[62];
	FFTPtr _fft;

	static std::vector<ElementPtr> elements;

	void setSymbol(std::string symbol)
	{
		_symbol = symbol;
	}

	void setName(std::string name)
	{
		_name = name;
	}
};

#endif /* defined(__vagabond__Element__) */
