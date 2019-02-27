//
//  Vagabond2GL.cpp
//  VagabondViewer
//
//  Created by Helen Ginn on 03/11/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#include "Vagabond2GL.h"
#include "../libsrc/Options.h"
#include "../libsrc/Crystal.h"
#include "../libsrc/Bond.h"
#include "../libsrc/Atom.h"
#include "../libsrc/Element.h"

bool Vagabond2GL::isAcceptableAtom(Atom *atom)
{
	if (!atom)
	{
		return false;
	}

	if (!atom->getElement())
	{
		return false;
	}

	if (atom->getElectronCount() <= 1)
	{
		return false;
	}

	if (atom->getModel() && !atom->getModel()->isBond())
	{
		return false;
	}

	return true;
}

bool Vagabond2GL::shouldGetBonds()
{
	if (!_moleculeMap.size() || _shouldGetBonds)
	{
		return true;
	}

	_renders = 0;
	
	/* Have there been any changes? */

	OptionsPtr globalOptions = Options::getRuntimeOptions();

	for (int i = 0; i < globalOptions->crystalCount(); i++)
	{
		CrystalPtr crystal = globalOptions->getCrystal(i);

		for (int j = 0; j < crystal->moleculeCount(); j++)
		{
			MoleculePtr molecule = crystal->molecule(j);

			if (!_moleculeMap.count(molecule))
			{
				return true;
			}

			int expected = _moleculeMap[molecule];
			int existing = 0;

			for (int i = 0; i < molecule->atomCount(); i++)
			{
				AtomPtr minor = molecule->atom(i);

				if (!isAcceptableAtom(&*minor))
				{
					continue;
				}
				
				if (!acceptablePositions(minor))
				{
					continue;
				}

				existing++;
			}
			
			if (expected != existing)
			{
				return true;
			}
		}
	}

	return false;
}

void Vagabond2GL::setVertexColour(AtomPtr atom, Vertex *vertex)
{
	vertex->color[0] = 100. / 255.;
	vertex->color[1] = 100. / 255.;
	vertex->color[2] = 100. / 255.;
	vertex->color[3] = 1.0;
	
	std::string symbol = atom->getElementSymbol();

	if (symbol == "O")
	{
		vertex->color[0] = 1.0;
		vertex->color[1] = 0.0;
		vertex->color[2] = 0.0;
	}

	if (symbol == "N")
	{
		vertex->color[0] = 104. / 255.;
		vertex->color[1] = 139. / 255.;
		vertex->color[2] = 255. / 255.;
	}

	if (symbol == "S")
	{
		vertex->color[0] = 255. / 255.;
		vertex->color[1] = 255. / 255.;
		vertex->color[2] = 0. / 255.;
	}
}

void Vagabond2GL::findAtoms()
{
	OptionsPtr globalOptions = Options::getRuntimeOptions();

	clearVertices();
	_moleculeMap.clear();

	if (globalOptions)
	{
		for (int i = 0; i < globalOptions->crystalCount(); i++)
		{
			CrystalPtr crystal = globalOptions->getCrystal(i);

			for (int j = 0; j < crystal->moleculeCount(); j++)
			{
				MoleculePtr molecule = crystal->molecule(j);

				int expected = processMolecule(molecule);

				_moleculeMap[molecule] = expected;
			}
		}
	}
}

void Vagabond2GL::render()
{
	if (shouldGetBonds() && !_pause)
	{
		findAtoms();
	}
	else if (!_pause)
	{
		updateAtoms();
	}

	rebindProgram();
}

AtomPtr Vagabond2GL::findAtomAtXY(double x, double y, double *z)
{
	vec3 target = make_vec3(x, y, 0);
	AtomPtr chosen = AtomPtr();

	for (int i = 0; i < _pairList.size(); i++)
	{
		if (_pairList[i].min.expired())
		{
			continue;
		}

		AtomPtr atom = _pairList[i].min.lock();
		vec3 pos = atom->getAbsolutePosition();
		
		double last = 1;
		vec3 model = mat4x4_mult_vec3(modelMat, pos, &last);
		vec3 proj = mat4x4_mult_vec3(projMat, model, &last);
		
		vec3_mult(&proj, 1 / last);

		if (proj.x < -1 || proj.x > 1)
		{
			continue;
		}

		if (proj.y < -1 || proj.y > 1)
		{
			continue;
		}
		
		if (model.z > 0)
		{
			continue;
		}

		vec3 diff = vec3_subtract_vec3(proj, target);
		
		if (fabs(diff.x) < 0.02 && fabs(diff.y) < 0.02)
		{
			if (model.z > *z)
			{
				*z = model.z;
				chosen = atom;
			}
		}
	}
	
	return chosen;
}

void Vagabond2GL::bindTextures()
{
	int num = 1;
	_textures.resize(num);

	glGenTextures(num, &_textures[0]);
	glBindTexture(GL_TEXTURE_2D, _textures[0]);
	checkErrors();
}
