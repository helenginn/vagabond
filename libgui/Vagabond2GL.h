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
#include "GLObject.h"
#include "shared_ptrs.h"
#include <map>
#include "Bond.h"

typedef std::map<AtomPtr, std::pair<int, int> > AtomMap;
typedef std::map<MoleculePtr, int> GLMoleculeMap;

class Vagabond2GL : public GLObject
{
public:
	Vagabond2GL()
	{
		_renders = 0;
	}
	void findAtoms();

	virtual void render();

private:
	int processMolecule(MoleculePtr molecule);
	vec3 _centroid;
	void updateAtoms();
	bool shouldGetBonds();

	void getPositions(AtomPtr atom, std::vector<vec3> *min,
					  std::vector<vec3> *maj);
	AtomMap _atomMap;
	GLMoleculeMap _moleculeMap;
	int _renders;
};

#endif /* Vagabond2GL_hpp */
