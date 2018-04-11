// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#ifndef __vagabond__VScript__
#define __vagabond__VScript__

#include <string>
#include <vector>

#include "shared_ptrs.h"

typedef enum
{
	VErrorInvalidThingType,
	VErrorThingRedeclaration,
	VErrorImmutableObject,
	VErrorLeftThingNotFound,
	VErrorReachedEOF,
	VErrorExpectedSemicolon,
	VErrorExpectedComma,
	VErrorExpectedBracket,
	VErrorGetterDoesNotExist,
	VErrorCounterDoesNotExist,
	VErrorCounterInappropriate,
	VErrorSimpleTypeFunctionCall,
	VErrorMissingImplementation,
	VErrorAssignmentOfVoid,
	VErrorOperationOnVoid,
	VErrorTypeMismatch,
	VErrorExpectedEquals,
	VErrorExpectedOperator,
	VErrorBeyondArrayBounds,
	VErrorMissingParameter,
	VErrorInappropriateParameter,
	VErrorInappropriateOperation,
	VErrorInappropriateScopeEnd,
} VScriptError;

typedef enum
{
	VCompEqual,
	VCompLessThan,
	VCompGreaterThan,
	VCompUnassigned,
} VScriptComparison;

/**
 * \class VScript
 * \brief Loads and executes VScript strings.
 */

class VScript
{
public:
	VScript();	
	
	void loadScript(std::string script);
	void execute();
private:
	std::string _script;	
	std::string _scriptCopy;	
	char *_char;

	void repairScript();
	bool parse();

	size_t charactersIn(char *pos);
	LeftThingPtr getLeftThing(std::string name);
	
	/* Decides if the future thing is a Thing. Send in a temporary
	* 	character position instead of the master. */
	bool isBetterThing(char *tmp);

	ThingPtr getThing(char **pos, ThingPtr left = ThingPtr(), 
	                  bool defRight = false);
	ThingPtr processRest(char **_char, ThingPtr rightThing);
	LeftThingPtr createLeftThing(std::string type, std::string name);

	ThingPtr getStringThing(char **pos);
	ThingPtr getNumberThing(char **pos);

	void skipNextScope(char **_char);
	void executeNextScope(char **_char);

	void makeNewScope();
	void loseScope();
	void reportLine();

	bool evaluateCondition(char **_char);
	void validate(char *pos);
	void incrementAndValidate(char **pos);
	std::string getNextWord(char **white, char limit = '\0');
	void handleError(VScriptError error);

	VScopePtr currentScope()
	{
		return _scopes[_scopes.size() - 1];
	}
	
	/* Will act like a stack */
	std::vector<VScopePtr> _scopes;
};

#endif
