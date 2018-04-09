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
	VErrorLeftThingNotFound,
	VErrorReachedEOF,
	VErrorExpectedSemicolon,
	VErrorGetterDoesNotExist,
	VErrorCounterDoesNotExist,
	VErrorCounterInappropriate,
	VErrorSimpleTypeFunctionCall,
	VErrorMissingImplementation,
	VErrorAssignmentOfVoid,
	VErrorTypeMismatch,
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
	void parse();

	size_t charactersIn(char *pos);
	LeftThingPtr getLeftThing(std::string name);
	
	/* Decides if the future thing is a Thing. Send in a temporary
	* 	character position instead of the master. */
	bool isThing(char *white, char *tmp);
	ThingPtr getThing(char **pos);
	LeftThingPtr createLeftThing(std::string type, std::string name);

	void makeNewScope();
	void loseScope();

	void validate(char *pos);
	void incrementAndValidate(char **pos);
	char *nextWhiteValidate(char *pos);
	void wrapNextWord(char **pos, char **white);

	VScopePtr currentScope()
	{
		return _scopes[_scopes.size() - 1];
	}
	
	/* Will act like a stack */
	std::vector<VScopePtr> _scopes;
};
