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

typedef std::map<AtomWkr, std::pair<int, int> > AtomMap;
typedef std::map<MoleculePtr, int> GLMoleculeMap;

class Vagabond2GL : public GLObject
{
public:
	Vagabond2GL()
	{
		_renders = 0;
		_lastEnsembleCount = 0;
		_shouldGetBonds = true;
		
		_pause = false;
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
	
	virtual void pause(bool on)
	{
		_pause = on;
		_shouldGetBonds = true;
	}
protected:
	virtual void bindTextures() = 0;
	virtual void updateAtoms() = 0;
	virtual void getPositions(AtomPtr atom, 
	                          std::vector<vec3> *min,
	                          std::vector<vec3> *maj) = 0;
	virtual int processMolecule(MoleculePtr molecule) = 0;

	void setVertexColour(AtomPtr atom, Vertex *vertex);
	int _lastEnsembleCount;
	bool _shouldGetBonds;

	AtomMap _atomMap;
	GLMoleculeMap _moleculeMap;
private:
	vec3 _centroid;
	bool shouldGetBonds();


	int _renders;
	bool _enabled;
	bool _pause;

};

#endif /* Vagabond2GL_hpp */
