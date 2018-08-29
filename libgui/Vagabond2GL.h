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
#include "../libsrc/shared_ptrs.h"
#include <map>
#include "../libsrc/Bond.h"

typedef std::map<AtomPtr, std::pair<int, int> > AtomMap;
typedef std::map<MoleculePtr, int> GLMoleculeMap;

class Vagabond2GL : public GLObject
{
public:
	Vagabond2GL(int average = false)
	{
		_renders = 0;
		_average = average;
		setupAverage();
		
		_enabled = true;
	}

	void findAtoms();

	virtual void render();
	
	bool isEnabled()
	{
		return _enabled;
	}
	
	void setEnabled(bool enabled)
	{
		_enabled = enabled;
	}
protected:
	virtual void bindTextures();
private:
	int processMolecule(MoleculePtr molecule);
	vec3 _centroid;
	void updateAtoms();
	bool shouldGetBonds();
	void setVertexColour(AtomPtr atom, Vertex *vertex);

	void getPositions(AtomPtr atom, std::vector<vec3> *min,
					  std::vector<vec3> *maj);
	AtomMap _atomMap;
	GLMoleculeMap _moleculeMap;
	void setupAverage();
	int _renders;
	int _average;
	bool _enabled;

};

#endif /* Vagabond2GL_hpp */
