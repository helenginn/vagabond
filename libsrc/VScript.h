// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

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
	VErrorGetterDoesNotExist,
	VErrorCounterDoesNotExist,
	VErrorCounterInappropriate,
	VErrorSimpleTypeFunctionCall,
	VErrorMissingImplementation,
	VErrorAssignmentOfVoid,
	VErrorOperationOnVoid,
	VErrorTypeMismatch,
	VErrorExpectedEquals,
	VErrorBeyondArrayBounds,
	VErrorMissingParameter,
	VErrorInappropriateOperation,
} VScriptError;

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
	bool isThing(char *tmp, char *white, bool nowind = false);
	ThingPtr getThing(char **pos, bool nowind = false, ThingPtr left = ThingPtr());
	ThingPtr processRest(char **_char, ThingPtr rightThing);
	LeftThingPtr createLeftThing(std::string type, std::string name);

	ThingPtr getStringThing(char **pos);
	ThingPtr getNumberThing(char **pos);

	void makeNewScope();
	void loseScope();
	void reportLine();

	void validate(char *pos);
	void incrementAndValidate(char **pos);
	char *nextWhiteValidate(char *pos);
	void wrapNextWord(char **pos, char **white);
	void handleError(VScriptError error);

	VScopePtr currentScope()
	{
		return _scopes[_scopes.size() - 1];
	}
	
	/* Will act like a stack */
	std::vector<VScopePtr> _scopes;
};
