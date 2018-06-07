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
#include "Shouter.h"

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
	if (_scopes.size() <= 1)
	{
		throw VErrorInappropriateScopeEnd;
	}
	
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

/* Will put the 'white' character to the one after the word ending */
std::string VScript::getNextWord(char **white, char limit)
{
	char *pos = *white;
	incrementAndValidate(&pos);
	if (limit == '\0')
	{
		*white = strchrwhite(pos);
	}
	else
	{
		*white = strchr(pos, limit);	
	}

	validate(*white);

	char tmp = **white;
	**white = '\0';
	std::string word = pos;
	**white = tmp;
	
	/* if limit is not a white character, it cannot be automatically
	* incremented, so we must push it forwards one */
	
	if (limit != '\0')
	{
		(*white)++;
	}

	return word;
}

inline bool looksLikeNumber(char *tmp)
{
	if (*tmp >= '0' && *tmp <= '9')
	{
		return true;
	}
	else if (*tmp == '.' || *tmp == '-')
	{
		return true;
	}
	
	return false;
}

inline bool looksLikeString(char *tmp)
{
	if (*tmp == '"')
	{
		return true;
	}

	return false;
}

bool VScript::isBetterThing(char *tmp)
{
	std::string firstWord = getNextWord(&tmp);
	
	if (looksLikeNumber(&firstWord[0]))
	{
		return false;
	}
	
	if (looksLikeString(&firstWord[0]))
	{
		return false;
	}
	
	int pos = 0;
	int brackstack = 0;
	
	while (pos < firstWord.length())
	{
		if (firstWord[pos] == '(') brackstack++;
		if (firstWord[pos] == ')') brackstack--;
		
		/* If we are completing a previous bracket stack,
		 * nope. */
		if (brackstack < 0) return false;
		
		if (firstWord[pos] == '.')
		{
			return true;
		}
		
		pos++;	
	}
	
	return false;
}

/* Condition reads like: (variable > 0) with brackets */
bool VScript::evaluateCondition(char **_char)
{
	(*_char)++;
	incrementAndValidate(_char);

	if (**_char != '(')
	{
		throw VErrorExpectedBracket;
	}
	
	(*_char)++;
	incrementAndValidate(_char);
	ThingPtr firstSide = getThing(_char);
	
	std::string op = getNextWord(_char);
	
	if (op != "==" && op != ">" && op != "<")
	{
		throw VErrorExpectedOperator;
	}
	
	VScriptComparison comp;

	if (op == "==")
	{
		comp = VCompEqual;
	}
	else if (op == ">")
	{
		comp = VCompGreaterThan;
	}
	else if (op == "<")
	{
		comp = VCompLessThan;
	}
	
	ThingPtr secondSide = getThing(_char);

	if (**_char != ')')
	{
		throw VErrorExpectedBracket;
	}
	
	(*_char)++;
	incrementAndValidate(_char);
	
	VScriptComparison actual = firstSide->compareToThing(secondSide);
	
	return (actual == comp);
}

void VScript::skipNextScope(char **_char)
{
	if (**_char != '{')
	{
		throw VErrorExpectedBracket;
	}
	
	int brackets = 1;
	
	while (brackets > 0)
	{
		(*_char)++;
		validate(*_char);
		
		if (**_char == '{')
		{
			brackets++;
		}
		else if (**_char == '}')
		{
			brackets--;
		}
	}
	
	(*_char)++;
	validate(*_char);
}

void VScript::executeNextScope(char **_char)
{
	if (**_char != '{')
	{
		throw VErrorExpectedBracket;
	}
	
	(*_char)++;
	incrementAndValidate(_char);

	makeNewScope();
	
	bool more = true;

	while (more)
	{
		more = parse();
	}
	
	(*_char)++;
	validate(*_char);
}

bool VScript::parse()
{
	std::string word;

	char *first = _char;
	try
	{
		word = getNextWord(&_char);
	}
	catch (const VScriptError &error)
	{
		if (error == VErrorReachedEOF)
		{
			return false;
		}
	}
	
	/* check reserved keywords */
	if (word == "}")
	{
		loseScope();
		return false;
	}
	if (word == "if")
	{
		bool proceed = evaluateCondition(&_char);
		
		if (!proceed)
		{
			skipNextScope(&_char);
		}
		else
		{
			executeNextScope(&_char);
		}

		return true;
	}
	else if (word == "while")
	{
		char *first = _char;

		bool proceed = evaluateCondition(&_char);
		
		if (!proceed)
		{
			skipNextScope(&_char);
		}
		else
		{	
			while (proceed)
			{
				executeNextScope(&_char);
				_char = first;
				proceed = evaluateCondition(&_char);
			}

			skipNextScope(&_char);
		}

		return true;
	}

	/* must be expecting a left thing ... or up to two words. */

	/* we make a temporary char because we need to rewind if it's a
	* Thing. */

	bool right = isBetterThing(first);


	if (right)
	{
		_char = first;
		getThing(&_char);
		goto check_semicolon;
	}
	else
	{
		std::string firstWord = word;//std::string(_char);

		/* It wasn't suitable, so we carry on with tmp as default and
		* 	find out the next word. */

		std::string secondWord = getNextWord(&_char);
		
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
			
			std::string thirdWord = getNextWord(&_char);

			/* Confirm an equals sign in the middle */
			if (thirdWord != "=")
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
	
	(_char)++;
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
	
	for (size_t i = 0; i < errorPlace; i++)
	{
		std::cout << " ";
	}
	
	std::cout << "^ - here" << std::endl;
	std::cout << std::endl;
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
		
		case VErrorExpectedBracket:
		std::cout << "Expected a bracket" << std::endl;
		break;
		
		case VErrorExpectedEquals:
		std::cout << "Expected a = sign" << std::endl;
		break;
		
		case VErrorExpectedSemicolon:
		std::cout << "Expecting semicolon";
		break;
		
		case VErrorExpectedOperator:
		std::cout << "Expected comparison operator (==, < or >)";
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
		
		case VErrorInappropriateScopeEnd:
		std::cout << "Inappropriate end of scope";
		break;
		
		case VErrorExpectedComma:
		std::cout << "Expecting comma";
		break;
		
		case VErrorInappropriateParameter:
		std::cout << "Inappropriate parameter for function";
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
	incrementAndValidate(pos);
	char *first = *pos;
	bool isDouble = false;
	
	while (looksLikeNumber(*pos))
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
	incrementAndValidate(pos);
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
	
	char tmp = **pos;
	**pos = '\0';
	first++;

	std::string aString = first;
	**pos = tmp;

	ThingPtr thing = ThingPtr(new Thing());
	thing->setThingType(ThingString);
	thing->setStringValue(aString);
	
	(*pos)++;

	return thing;
}

/* Could be a function too, in which case may return ThingPtr(). */
ThingPtr VScript::getThing(char **pos, ThingPtr thing, bool defRight)
{
	char *white = *pos;

	/* Something is "right" if it involves a function on a LeftThing.
	* It is always false if "thing" is an object, as this will have
	* occurred at some point up the stack already. */
	bool right =  isBetterThing(*pos) || defRight;

	if (!right)
	{
		char *orig = *pos;
		std::string word = getNextWord(pos);

		/* Start with reserved things... like numbers, strings */
		if (looksLikeNumber(&word[0]))
		{
			*pos = orig;
			ThingPtr thing = getNumberThing(pos);
			incrementAndValidate(pos);
			thing = processRest(pos, thing);
			return thing;
		}
		else if (word[0] == '\"')
		{
			*pos = orig;
			ThingPtr thing = getStringThing(pos);
			incrementAndValidate(pos);
			thing = processRest(pos, thing);
			return thing;
		}

		/* Last option is it's a LeftThing masquerading as a Thing. */

		bool hasEndRubbish = false;
		bool done = false;

		while (!done)
		{
			char finalchar = word[word.length() - 1];
			if (finalchar == ';' || finalchar == ')' || 
			    finalchar == ',')
			{
				hasEndRubbish = true;
				word = word.substr(0, word.length() - 1);
				(*pos)--;
			}
			else
			{
				break;
			}
		}

		if (!hasEndRubbish)
		{
			incrementAndValidate(pos);
		}

		/* Place cursor at the semicolon when we return */
		LeftThingPtr left = getLeftThing(word);
		ThingPtr right = ToThingPtr(left);
		right = processRest(pos, right);
		return right;
	}
	
	/* Now it's a function being performed on a thing. */
	
	ThingPtr left = thing;
	if (!thing)
	{
		/* We need to get the thing from the right of the left dot */
		std::string word = getNextWord(pos, '.');
		left = getLeftThing(word);
	}
	
	std::string function = getNextWord(pos, '(');

	std::vector<ThingPtr> things;

	while (**pos != ')')
	{
		incrementAndValidate(pos);
		ThingPtr thing = getThing(pos);
		things.push_back(thing);
		
		if (**pos != ',' && **pos != ')')
		{
			throw VErrorExpectedComma;	
		}
		else if (**pos == ',')
		{
			(*pos)++;
		}
	}

	(*pos)++;
	incrementAndValidate(pos);	
	
	ThingPtr rightThing = left->dealWithFunction(function, things);	
	rightThing = processRest(pos, rightThing);
	return rightThing;
}

ThingPtr VScript::processRest(char **pos, ThingPtr rightThing)
{
	if (**pos == ';')
	{
		return rightThing;
	}
	else if (**pos == '.')
	{
		(*pos)++;
		rightThing = getThing(pos, rightThing, true);
	}
	else if (**pos == '+')
	{
		if (!rightThing)
		{
			throw VErrorOperationOnVoid;
		}
		(*pos)++;
		incrementAndValidate(pos);
		
		ThingPtr secondThing = getThing(pos);
		rightThing->addThing(secondThing);
		
		return rightThing;
	}
	
	return rightThing;
}

LeftThingPtr VScript::getLeftThing(std::string name)
{
	LeftThingPtr thing;

	std::string sanitised;
	for (int i = 0; i < name.length(); i++)
	{
		if ((name[i] >= 'A' && name[i] <= 'Z') || (name[i] == '_') ||
		    (name[i] >= 'a' && name[i] <= 'z') ||
			(name[i] >= '0' && name[i] <= '9' && i > 0))
		{
			sanitised.push_back(name[i]);	
		}
		else
		{
			break;
		}
	}

	/** Prioritise latest scopes first */
	for (int i = _scopes.size() - 1; i >= 0; i--)
	{
		VScopePtr scope = _scopes[i];
		
		thing = scope->findThing(sanitised);
		
		if (thing) break;
	}
	
	if (!thing)
	{
		warn_user("Cannot find " + sanitised);
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





