// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#include "VScript.h"
#include "VScope.h"
#include "charmanip.h"
#include "LeftThing.h"
#include "Thing.h"
#include "Options.h"

VScript::VScript()
{
	makeNewScope();

	LeftThingPtr left = createLeftThing("Crystal", "main_crystal");
	CrystalPtr crystal = Options::getRuntimeOptions()->getActiveCrystal();
	ParserPtr parser = ToParserPtr(crystal);
	
	left->setParserValue(parser);
}

void VScript::makeNewScope()
{
	/* Start in a scope. */
	VScopePtr scope = VScopePtr(new VScope());
	_scopes.push_back(scope);
}

void VScript::loseScope()
{
	_scopes.pop_back();
}

void VScript::loadScript(std::string script)
{
	_script = script;
}

void VScript::validate(char *pos)
{
	if (pos == NULL)
	{
		throw VErrorReachedEOF;
	}
	
	size_t size = charactersIn(pos);
	
	if (size > _scriptCopy.length())
	{
		throw VErrorReachedEOF;
	}
}	

void VScript::incrementAndValidate(char **pos)
{
	incrementIndent(pos);
	validate(*pos);
}

char *VScript::nextWhiteValidate(char *pos)
{
	char *white = strchrwhite(pos);
	validate(white);
	
	*white = '\0';
	return white;
}

/** pos will be placed at beginning of word, white will be placed at end of
 * word and set to NULL. Starting positions: *white at the end of the last
 * word. */
void VScript::wrapNextWord(char **pos, char **white)
{
	*pos = *white + 1;
	incrementAndValidate(&*pos);
	*white = nextWhiteValidate(*pos);
}

bool VScript::isThing(char *tmp, char *white, bool nowind)
{
	if (!nowind)
	{
		wrapNextWord(&tmp, &white);
	}

	std::string firstWord = std::string(tmp);

	if (!nowind)
	{
		*white = ' ';
	}
	
	if (firstWord[0] >= '0' && firstWord[0] <= '9')
	{
		return false;
	}
	
	return (firstWord.find('.') != std::string::npos);
}

bool VScript::parse()
{
	char *white;
	try
	{
		incrementAndValidate(&_char);
		white = nextWhiteValidate(_char);
	}
	catch (const VScriptError &error)
	{
		if (error == VErrorReachedEOF)
		{
			return false;
		}
	}
	
	/* check reserved keywords */
	if (_char == "if")
	{
		throw VErrorMissingImplementation;
	}
	else if (_char == "while")
	{
		throw VErrorMissingImplementation;
	}

	/* must be expecting a left thing ... or up to two words. */

	/* we make a temporary char because we need to rewind if it's a
	* Thing. */
	
	bool right = isThing(_char, white, true);

	if (right)
	{
		getThing(&_char, true);
		/* Ignore output */
		goto check_semicolon;
	}
	else
	{
		*white = '\0';
		std::string firstWord = std::string(_char);

		/* It wasn't suitable, so we carry on with tmp as default and
		* 	find out the next word. */

		wrapNextWord(&_char, &white);
		std::string secondWord = std::string(_char);
		
		/* Confirm an equals sign in the middle */

		if (secondWord == "=")
		{
			/* this is an existing LeftThing which is being assigned. */
			/* Get a right thing. */
			ThingPtr right = getThing(&_char);
			
			if (!right)
			{
				throw VErrorAssignmentOfVoid;	
			}
			
			LeftThingPtr left = getLeftThing(firstWord);
			
			left->setThing(right);
		}
		else
		{
			/* This might be the declaration of a LeftThing. */
			LeftThingPtr thing = createLeftThing(firstWord, secondWord);
			
			wrapNextWord(&_char, &white);
			/* Confirm an equals sign in the middle */
			std::string secondWord = std::string(_char);

			if (secondWord != "=")
			{
				throw VErrorExpectedEquals;
			}

			/* Get a right thing. */
			ThingPtr right = getThing(&_char);
			
			if (!right)
			{
				throw VErrorAssignmentOfVoid;
			}
			
			thing->setThing(right);
		}
	}
	
	/* Validate that we ended on a ; */

check_semicolon:	
	if (*_char != ';')
	{
		throw VErrorExpectedSemicolon;
	}
	
	_char++;
	return true;
}

void VScript::reportLine()
{
	size_t placeForward = charactersIn(_char);
	size_t placeBackward = charactersIn(_char);
	
	while (_scriptCopy[placeBackward] != '\n')
	{
		placeBackward--;
		if (placeBackward == 0) break;
	}
	
	if (placeBackward > 0)
	{	
		placeBackward++;
	}

	while (_scriptCopy[placeForward] != '\n')
	{
		placeForward++;
		if (placeForward == _scriptCopy.size() - 1) break;
	}

	size_t length = placeForward - placeBackward;
	int errorPlace = charactersIn(_char) - placeBackward;
	
	std::string line = _scriptCopy.substr(placeBackward, length);
	
	std::cout << "\t\t..." << std::endl;
	std::cout << "\t\t" << line << std::endl;

	if (errorPlace > length || errorPlace < 0)
	{
		return;
	}
	
	std::cout << "\t\t";
	
	for (size_t i = 0; i < errorPlace - 1; i++)
	{
		std::cout << " ";
	}
	
	std::cout << "^ - here" << std::endl;
	std::cout << std::endl;
	
	std::cout << _script << std::endl;
}

void VScript::handleError(VScriptError error)
{
	std::cout << std::endl;
	std::cout << "\t\t OH NO!! ****** I'm an error, handle me ;)" << std::endl;
	std::cout << "\t\t";

	switch (error)
	{
		case VErrorInvalidThingType:
		std::cout << "Invalid type specified";
		break;
		
		case VErrorThingRedeclaration:
		std::cout << "Attempting to redeclare pre-existing variable";
		break;
	
		case VErrorLeftThingNotFound:
		std::cout << "Using undefined variable";
		break;
		
		case VErrorReachedEOF:
		std::cout << "Unexpectedly reached end of file";
		break;
		
		case VErrorExpectedEquals:
		std::cout << "Expected a = sign." << std::endl;
		break;
		
		case VErrorExpectedSemicolon:
		std::cout << "Expecting semicolon";
		break;
		
		case VErrorGetterDoesNotExist:
		std::cout << "Object does not have this getter function";
		break;
	
		case VErrorCounterDoesNotExist:
		std::cout << "Object does not have this counter function";
		break;
	
		case VErrorCounterInappropriate:
		std::cout << "This property cannot be counted";
		break;
		
		case VErrorSimpleTypeFunctionCall:
		std::cout << "Attempting to call inappropriate function"\
		" on simple type (int, double etc.)";
		break;
		
		case VErrorBeyondArrayBounds:
		std::cout << "Beyond bounds of array." << std::endl;
		break;
		
		case VErrorOperationOnVoid:
		std::cout << "Attempting operation, but left hand object is void";
		break;
		
		case VErrorInappropriateOperation:
		std::cout << "Inappropriate operation on object types";
		break;
		
		case VErrorAssignmentOfVoid:
		std::cout << "Function returns nothing, attempting to assign to variable";
		break;
		
		case VErrorTypeMismatch:
		std::cout << "Return of function and variable have different types.";
		break;
		
		case VErrorMissingParameter:
		std::cout << "Missing parameter to function.";
		break;
	
		case VErrorMissingImplementation:
		std::cout << "Helen hasn't written code to support this yet";
		break;
	
		default:
		break;	
	}

	reportLine();
	std::cout << std::endl << std::endl;
}

void VScript::execute()
{
	std::cout << "***************************" << std::endl;
	std::cout << "**** Executing VScript ****" << std::endl;
	std::cout << "***************************" << std::endl;
	_script.push_back('\n');
	_scriptCopy = _script;
	_char = &_script[0];
	
	bool more = true;

	try
	{
		while (more)
		{
			more = parse();
		}
	}
	catch (VScriptError &error)
	{
		handleError(error);
	}
	
	std::cout << "***************************" << std::endl;
	std::cout << "****       End.         ****" << std::endl;
	std::cout << "***************************" << std::endl;
	
	repairScript();
}

size_t VScript::charactersIn(char *pos)
{
	size_t difference = pos - &_script[0];
	return difference;
}

void VScript::repairScript()
{
	_script = _scriptCopy;
}

ThingPtr VScript::getNumberThing(char **pos)
{
	char *first = *pos;
	bool isDouble = false;
	
	while ((**pos >= '0' && **pos <= '9') || **pos == '.')
	{
		if (**pos == '.')
		{
			isDouble = true;
		}
		
		(*pos)++;
		validate(*pos);
	}
	
	/* Get word without disturbing much */
	char store = (**pos);
	**pos = '\0';
	std::string word = std::string(first);
	**pos = store;
	
	ThingPtr thing = ThingPtr(new Thing());
	
	if (isDouble)
	{
		double dval = atof(word.c_str());
		thing->setThingType(ThingDouble);
		thing->setDoubleValue(dval);
	}
	else
	{
		int ival = atoi(word.c_str());
		thing->setThingType(ThingInt);
		thing->setIntValue(ival);
	}
	
	return thing;
}

ThingPtr VScript::getStringThing(char **pos)
{
	bool freePass = false;
	char *first = *pos;
	(*pos)++;
	
	while (**pos != '"' || freePass)
	{
		(*pos)++;

		validate(*pos);
		
		if (**pos == '\\' && freePass == false)
		{
			freePass = true;
		}
		else
		{
			freePass = false;
		}
	}	
	
	**pos = '\0';
	*first++;

	std::string aString = first;
	ThingPtr thing = ThingPtr(new Thing());
	thing->setThingType(ThingString);
	thing->setStringValue(aString);
	
	(*pos)++;

	return thing;
}

/* Could be a function too, in which case may return ThingPtr(). */
ThingPtr VScript::getThing(char **pos, bool nowind, ThingPtr thing)
{
	char *white = *pos;

	/* Something is "right" if it involves a function on a LeftThing.
 	 * It is always false if "thing" is an object, as this will have
	 * occurred at some point up the stack already. */
	bool right = (thing != ThingPtr());

	/* If Thing is not set, however, we should check to see if it is
	 * a right thing. */
	if (!thing)
	{
		char *tmp = *pos;
		right = isThing(tmp, white, nowind);
	}

	if (!nowind)
	{
		white = *pos;
		wrapNextWord(pos, &white);
	}

	/* Start with reserved things... like numbers, strings */
	if (**pos >= '0' && **pos <= '9')
	{
		ThingPtr thing = getNumberThing(pos);
		thing = processRest(pos, thing);
		return thing;
	}
	else if (**pos == '\"')
	{
		/* Restore the white character again */
		*white = ' ';

		ThingPtr thing = getStringThing(pos);
		thing = processRest(pos, thing);
		return thing;
	}

	/* If we are here but it's not a function, it must be a leftThing */
	if (!right)
	{
		/* Last option is it's a LeftThing masquerading as a Thing. */
		/* n.b. we should have a trailing semicolon. */

		std::string word = std::string(*pos);
		*white = ' ';

		bool hasSemicolon = false;

		if (word[word.length() - 1] == ';')
		{
			hasSemicolon = true;
		}

		if (hasSemicolon)
		{
			word.pop_back();
		}

		/* Place cursor at the semicolon when we return */
		*pos += word.length();

		LeftThingPtr left = getLeftThing(word);
		ThingPtr right = ToThingPtr(left);
		right = processRest(pos, right);
		return right;
	}
	
	/* Now it's a function being performed on a thing. */
	
	char *dot = strchr(*pos, '.');
	validate(dot);
	*dot = '\0';
	
	ThingPtr left = thing;
	
	if (!thing)
	{
		/* word should now contain the left thing */
		std::string word = std::string(*pos);
		left = getLeftThing(word);
	}
	
	*pos = dot + 1;
	validate(*pos);
	
	char *leftBracket = strchr(*pos, '(');
	validate(leftBracket);
	*leftBracket = '\0';
	std::string function = std::string(*pos);
	
	leftBracket++;
	char *rightBracket = strchr(leftBracket, ')');
	validate(rightBracket);
	*rightBracket = '\0';
	std::string contents = std::string(leftBracket);
	
	*pos = rightBracket + 1;

	ThingPtr rightThing = left->dealWithFunction(function, contents);	
	rightThing = processRest(pos, rightThing);
	return rightThing;
}

ThingPtr VScript::processRest(char **_char, ThingPtr rightThing)
{
	incrementAndValidate(_char);
	
	if (**_char == ';')
	{
		return rightThing;
	}
	else if (**_char == '.')
	{
		rightThing = getThing(_char, true, rightThing);
	}
	else if (**_char == '+')
	{
		if (!rightThing)
		{
			throw VErrorOperationOnVoid;
		}
//		(*_char)++;
//		incrementAndValidate(_char);
		
		ThingPtr secondThing = getThing(_char);
		
		rightThing->addThing(secondThing);
	}
	
	return rightThing;
}

LeftThingPtr VScript::getLeftThing(std::string name)
{
	LeftThingPtr thing;

	/** Prioritise latest scopes first */
	for (int i = _scopes.size() - 1; i >= 0; i--)
	{
		VScopePtr scope = _scopes[i];
		
		thing = scope->findThing(name);
		
		if (thing) break;
	}
	
	if (!thing)
	{
		throw VErrorLeftThingNotFound;
	}
	
	return thing;
}

LeftThingPtr VScript::createLeftThing(std::string type, std::string name)
{
	LeftThingPtr left = LeftThingPtr(new LeftThing());
	left->setThingType(type);
	left->setThingName(name);

	currentScope()->addLeftThing(left);
	return left;
}





