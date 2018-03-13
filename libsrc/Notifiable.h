//
//  Notifiable.h
//  Notifiable
//
//  Created by Helen Ginn on 21/01/2018
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#ifndef __Vagabond__Notifiable__
#define __Vagabond__Notifiable__

#include <string>
#include "AtomGroup.h"

typedef enum
{
	InstructionTypeNone,
	InstructionTypeOpenPDB,
	InstructionTypeOpenMTZ,
	InstructionTypeSuperimpose,
	InstructionTypeRefinePositions,
	InstructionTypeRefineFlexibility,
	InstructionTypeRefineDensity,
	InstructionTypeSidechainsToEnd,
	InstructionTypeRefineToEnd,
	InstructionTypeChangeBMult,
	InstructionTypeRecalculateFFT,
	InstructionTypeSetOutputDir,
	InstructionTypeSetObjectValue,
	InstructionTypeGetObjectValue,
	InstructionTypeFitWholeMoleculeTranslation,
	InstructionTypeFitWholeMoleculeRotation,
} InstructionType;

class Notifiable
{
public:
	Notifiable()
	{
		_value = 0;
		_result = 0;
		_setter = NULL;
		_getter = NULL;
		_object = NULL;
		_enabled = false;
		_atomGroup = AtomGroupPtr();
	}

	virtual void enable()
	{
		_enabled = true;
	}

	virtual void disable()
	{
		_enabled = false;
	}

	virtual bool isRunningSomething() = 0;

	virtual void setMessage(std::string message)
	{
		_message = message;
	}

	void setObject(void *object)
	{
		_object = object;
	}

	void setSetter(Setter setter, double value)
	{
		_value = value;
		_setter = setter;
	}    

	void setGetter(Getter getter)
	{
		_getter = getter;
	}

	void performObjectSet()
	{
		disable();
		(*_setter)(_object, _value);

		if (_atomGroup)
		{
			_atomGroup->refreshPositions(false);
		} 

		enable();
	}

	void performObjectGet()
	{
		disable();
		_result = (*_getter)(_object);

		if (_atomGroup)
		{
			_atomGroup->refreshPositions(true);
		}

		enable();
	}

	void setInstruction(InstructionType type)
	{
		_instructionType = type;
		wakeup();
	}

	virtual void wakeup() = 0;

	void setRefreshGroup(AtomGroupPtr group)
	{
		_atomGroup = group;
	}
protected:
	InstructionType _instructionType;

	void *getObject()
	{
		return _object;
	}

private:
	std::string _message;
	bool _enabled;
	void *_object;
	Setter _setter;
	Getter _getter;
	double _value;
	double _result;
	AtomGroupPtr _atomGroup;
};

#endif
