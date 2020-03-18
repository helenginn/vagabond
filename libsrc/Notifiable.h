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
	InstructionTypeResetExplorer,
	InstructionTypeRefinePositions,
	InstructionTypeRefinePosToDensity,
	InstructionTypeRefineFlexibility,
	InstructionTypeRefineDensity,
	InstructionTypeSidechainsToEnd,
	InstructionTypeRefineToEnd,
	InstructionTypeModelPosToEnd,
	InstructionTypeRecalculateFFT,
	InstructionTypeOmitScan,
	InstructionTypeOpenInCoot,
	InstructionTypeSetObjectValue,
	InstructionTypeGetObjectValue,
	InstructionTypePreviousState,
	InstructionTypeRefineSidePos,
	InstructionTypeRigidBody,
	InstructionTypeAddPDBFile,
	InstructionTypeUpdateFromPDB,
	InstructionTypeFitTranslation,
	InstructionTypeResetMotion,
	InstructionTypeResetSides,
	InstructionTypeResetIntra,
	InstructionTypeFindDisulphides,
	InstructionTypeRefineIntramolecule,
	InstructionTypeSplitBond,
	InstructionTypeAdjustBFactor,
	InstructionTypeRefitBackbone,
	InstructionTypeChelate,
	InstructionTypeSponge,
	InstructionTypeCancel,
	InstructionTypeManualRefine,
	InstructionTypeWritePNG,
	InstructionTypeRefineIntramagic,
} InstructionType;

/**
 * \class Notifiable
 * \brief Allows interface to a GUI. This should be overloaded by a class in
 * the GUI in order to respond to changes from Vagabond during refinement.
 */

class Notifiable
{
public:
	Notifiable()
	{
		_value = 0;
		_cancel = false;
		_result = 0;
		_setter = NULL;
		_getter = NULL;
		_object = NULL;
		_enabled = false;
		_atomGroup = AtomGroupPtr();
	}
	
	virtual ~Notifiable()
	{

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
	
	void setValue(double value)
	{
		_value = value;
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
	
	bool isEnabled()
	{
		return _enabled;
	}

	void setRefreshGroup(AtomGroupPtr group)
	{
		_atomGroup = group;
	}
	
	virtual void focusOnPosition(vec3 pos)
	{

	}
	
	virtual void appendToLog(char *msg) = 0;
	virtual void pause(bool on) = 0;
	
	bool shouldCancel()
	{
		return _cancel;
	}
	
	void setDoneCancelling()
	{
		_cancel = false;
	}

	virtual void setRenderDensity() = 0;
	virtual void renderWarp() = 0;
protected:
	InstructionType _instructionType;

	void *getObject()
	{
		return _object;
	}
	
	double getValue()
	{
		return _value;
	}

	void setShouldCancel()
	{
		_cancel = true;
		std::cout << "Should cancel after next FFT." << std::endl;
		setMessage("Cancelling after next FFT.");
	}
private:
	std::string _message;
	bool _enabled;
	bool _cancel;
	void *_object;
	Setter _setter;
	Getter _getter;
	double _value;
	double _result;
	AtomGroupPtr _atomGroup;
};

#endif
