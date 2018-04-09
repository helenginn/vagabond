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
 * word and set to NULL. */
void VScript::wrapNextWord(char **pos, char **white)
{
	*pos = *white + 1;
	incrementAndValidate(&*pos);
	*white = nextWhiteValidate(*pos);
}

bool VScript::isThing(char *white, char *tmp)
{
	std::string firstWord = std::string(tmp);
	
	if (firstWord.find('.') != std::string::npos)
	{
		/* this is a LeftThing for which a function is being called.
		* Reset situation ... */
		*white = ' ';
		
		return true;
	}
	
	return false;
}

void VScript::parse()
{
	incrementAndValidate(&_char);

	char *white = nextWhiteValidate(_char);
	
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
	char *tmp = _char;
	wrapNextWord(&tmp, &white);
	
	bool right = isThing(white, tmp);

	if (right)
	{
		getThing(&_char);
		/* Ignore output */
		goto check_semicolon;
	}
	else
	{
		std::string firstWord = std::string(_char);

		/* It wasn't suitable, so we carry on with tmp as default and
		* 	find out the next word. */
		_char = tmp;

		wrapNextWord(&_char, &white);
		std::string secondWord = std::string(_char);

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
			
			/* Confirm an equals sign in the middle */
			throw VErrorMissingImplementation;

			/* Get a right thing. */
			getThing(&_char);
		}
	}
	
	/* Validate that we ended on a ; */

check_semicolon:	
	if (*_char != ';')
	{
		throw VErrorExpectedSemicolon;
	}
}

void VScript::execute()
{
	_scriptCopy = _script;
	_char = &_script[0];

	try
	{
		parse();
	}
	catch (const VScriptError)
	{
		std::cout << "I'm an error, handle me ;)" << std::endl;
	}
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

/* Could be a function too, in which case may return ThingPtr(). */
ThingPtr VScript::getThing(char **pos)
{
	/* Start with reserved things... like numbers, strings */
	char *tmp = *pos;
	char *white = *pos;
	wrapNextWord(&tmp, &white);
	
	if (**pos >= '0' && **pos <= '9')
	{
		/* We have an int or a double? */
		throw VErrorMissingImplementation;
	}
	else if (**pos == '\"')
	{
		/* We have a string? */
		throw VErrorMissingImplementation;
	}
	
	bool right = isThing(white, tmp);
	
	if (!right)
	{
		/* Last option is it's a LeftThing masquerading as a Thing. */
		/* n.b. we should have a trailing semicolon. */
		
		std::string word = std::string(*pos);
		if (word[word.length() - 1] != ';')
		{
			throw VErrorMissingImplementation;
		}
		
		/* Remove last semicolon */
		word.pop_back();
		
		/* Place cursor at the semicolon when we return */
		*pos += word.length();
		
		LeftThingPtr left = getLeftThing(word);
		ThingPtr rightEquiv = ToThingPtr(left);
		return rightEquiv;
	}
	
	/* Now it's a function being performed on a left thing. */
	
	char *dot = strchr(*pos, '.');
	validate(dot);
	*dot = '\0';
	
	/* word should now contain the left thing */
	std::string word = std::string(*pos);
	LeftThingPtr left = getLeftThing(word);
	
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

	ThingPtr rightThing = left->dealWithFunction(function, contents);	
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
}





