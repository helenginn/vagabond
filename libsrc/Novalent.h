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

#ifndef __vagabond__novalent__
#define __vagabond__novalent__

#include "libsrc/shared_ptrs.h"
#include "ExplicitModel.h"

class Novalent : public ExplicitModel
{
public:
	Novalent();
	Novalent(AtomPtr atom);
	
	virtual ~Novalent();

	virtual std::vector<BondSample> *getManyPositions(void *object = NULL);
	virtual void propagateChange(int depth = -1, bool refresh = false);
	virtual std::string shortDesc();
	virtual AnchorPtr getAnchor();

	virtual std::string getClassName()
	{
		return "Novalent";
	}
	
	virtual AtomPtr getAtom()
	{
		return _atom.lock();
	}
	
	virtual bool isDisabled()
	{
		return false;
	}
protected:
	void initialise();
	virtual std::string getParserIdentifier();
	virtual void addProperties();
	virtual void postParseTidy();
	virtual void linkReference(BaseParserPtr object, std::string category); 
	virtual void addObject(ParserPtr object, std::string category);

	vec3 _abs;
private:
	AnchorPtr _anch;
	AtomWkr _atom;
	double _b;
};

#endif
