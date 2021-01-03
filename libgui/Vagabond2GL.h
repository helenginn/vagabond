//
//  Vagabond2GL.hpp
//  VagabondViewer
//
//  Created by Helen Ginn on 03/11/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#ifndef Vagabond2GL_hpp
#define Vagabond2GL_hpp

#include <stdio.h>
#include "../subprojects/helen3d/libsrc/SlipObject.h"
#include "../libsrc/shared_ptrs.h"
#include <map>
#include "../libsrc/Bond.h"

typedef struct
{
	AtomWkr maj;
	AtomWkr min;
	int vNum;
	int size;
} Atom3D;

typedef std::map<MoleculePtr, int> GLMoleculeMap;

class Vagabond2GL : public SlipObject
{
public:
	Vagabond2GL();

	virtual void render(SlipGL *sender);
	
	bool isEnabled()
	{
		return _enabled;
	}
	
	void setEnabled(bool enabled)
	{
		_enabled = enabled;
	}
	
	virtual void pause(bool on)
	{
		_pause = on;
		_shouldGetBonds = true;
	}
	
	AtomPtr findAtomAtXY(double x, double y, double *z);
	
	virtual bool isVagabond2GL()
	{
		return true;
	}
	
	void setShouldGetBonds(bool should = true)
	{
		_shouldGetBonds = should;
	}
	
	void toggleColourByFlex()
	{

	}
protected:
	virtual void findAtoms();
	bool isAcceptableAtom(Atom *atom);

	virtual void extraUniforms();
	virtual void updateAtoms() = 0;
	void updateColour(AtomPtr atom, Helen3D::Vertex *vertex);
	virtual bool acceptablePositions(AtomPtr minAtom)
	{
		return true;
	};

	virtual bool getPositions(AtomPtr minAtom, AtomPtr majAtom, 
	                          std::vector<vec3> *min,
	                          std::vector<vec3> *maj) = 0;
	virtual int processMolecule(MoleculePtr molecule) = 0;

	void setVertexColour(AtomPtr atom, Helen3D::Vertex *vertex);
	int _lastEnsembleCount;
	bool _shouldGetBonds;
	bool _enabled;
	bool _colourByFlex;

	/* First pair: number of bonds assoc. with atom,
	 * Second pair: vertex num */
	std::vector<Atom3D> _pairList;
	GLMoleculeMap _moleculeMap;
private:
	vec3 _centroid;
	bool shouldGetBonds();

	int _renders;
	bool _pause;
};

#endif /* Vagabond2GL_hpp */
