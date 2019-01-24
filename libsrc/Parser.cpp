// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#include "Parser.h"
#include <iostream>
#include <iomanip>
#include "Crystal.h"
#include "Polymer.h"
#include "Twist.h"
#include "WaterNetwork.h"
#include "Atom.h"
#include "Bond.h"
#include "BondGroup.h"
#include "Monomer.h"
#include "Absolute.h"
#include "Anchor.h"
#include "Sidechain.h"
#include "Backbone.h"
#include "Shouter.h"
#include "charmanip.h"

ParserMap Parser::_allParsers;
ClassMap Parser::_allClasses;

Parser::Parser()
{
	_restored = false;
	_setup = false;
	_parent = NULL;
	setupKnownClasses();
}

void Parser::setup(bool isNew)
{
	if (_setup) return;

	/* If these objects have been generated during this session,
	 * not loaded from a file... */
	if (!isNew)
	{
		_identifier = getParserIdentifier(); 

		makePath();
		
		std::string path = getAbsolutePath();
		ParserPtr newPointer = shared_from_this();
		addToAllParsers(path, newPointer);
	}

	_className = getClassName();

	addProperties();
}

void Parser::makePath()
{
	std::string path;

	if (_parent)
	{
		path = _parent->getAbsolutePath();
	}
	else
	{
		path = "";
	}

	path += "/";
	path += _className + "::" + _identifier;
	_absolutePath = path;
}

void Parser::setParent(Parser *parent)
{
	_parent = parent;
	makePath();
}

void Parser::addStringProperty(std::string className, std::string *ptr)
{
	StringProperty property;
	property.ptrName = className;
	property.stringPtr = ptr;
	_stringProperties.push_back(property);
}

std::string *Parser::getStringProperty(std::string className)
{
	std::vector<StringProperty>::iterator it;
	
	for (int i = 0; i < _stringProperties.size(); i++)
	{
		if (_stringProperties[i].ptrName == className)
		{
			return _stringProperties[i].stringPtr;	
		}
	}

	return NULL;
}

void Parser::addDoubleProperty(std::string className, double *ptr)
{
	DoubleProperty property;
	property.ptrName = className;
	property.doublePtr = ptr;
	_doubleProperties.push_back(property);
}

double *Parser::getDoubleProperty(std::string className)
{
	std::vector<DoubleProperty>::iterator it;
	
	for (int i = 0; i < _doubleProperties.size(); i++)
	{
		if (_doubleProperties[i].ptrName == className)
		{
			return _doubleProperties[i].doublePtr;	
		}
	}
	
	return NULL;
}

void Parser::addIntProperty(std::string className, int *ptr)
{
	IntProperty property;
	property.ptrName = className;
	property.intPtr = ptr;
	_intProperties.push_back(property);
}

int *Parser::getIntProperty(std::string className)
{
	std::vector<IntProperty>::iterator it;
	
	for (int i = 0; i < _intProperties.size(); i++)
	{
		if (_intProperties[i].ptrName == className)
		{
			return _intProperties[i].intPtr;	
		}
	}

	return NULL;
}

void Parser::addMat3x3Property(std::string className, mat3x3 *ptr)
{
	Mat3x3Property property;
	property.ptrName = className;
	property.mat3x3Ptr = ptr;
	_mat3x3Properties.push_back(property);
}

void Parser::addVec3Property(std::string className, vec3 *ptr)
{
	Vec3Property property;
	property.ptrName = className;
	property.vec3Ptr = ptr;
	_vec3Properties.push_back(property);
}

void Parser::addMat3x3ArrayProperty(std::string className, std::vector<mat3x3> *ptr)
{
	Mat3x3ArrayProperty property;
	property.ptrName = className;
	property.mat3x3ArrayPtr = ptr;
	_mat3x3ArrayProperties.push_back(property);
}

void Parser::addVec3ArrayProperty(std::string className, std::vector<vec3> *ptr)
{
	Vec3ArrayProperty property;
	property.ptrName = className;
	property.vec3ArrayPtr = ptr;
	_vec3ArrayProperties.push_back(property);
}

void Parser::addBoolProperty(std::string className, bool *ptr)
{
	BoolProperty property;
	property.ptrName = className;
	property.boolPtr = ptr;
	_boolProperties.push_back(property);
}

void Parser::addChild(std::string category, ParserPtr child)
{
	if (!child) return;
	
	child->setParent(this);
	_parserList[category].push_back(child);
	child->setup();
}

ParserPtr Parser::getChild(std::string className, int num)
{
	if (!_parserList.count(className))
	{
		return ParserPtr();	
	}

	return _parserList[className][num];
}

void Parser::exposeFunction(std::string funcName, TwoDouble func)
{
	_twoDoubleList[funcName] = func;
}

void Parser::exposeFunction(std::string funcName, Getter func)
{
	_functionList[funcName] = func;
}

void Parser::exposeFunction(std::string funcName, Setter func)
{
	_setterList[funcName] = func;
}

int Parser::getChildCount(std::string className)
{
	if (!_parserList.count(className))
	{
		return 0;	
	}

	return _parserList[className].size();
}

void Parser::sanitise(std::string *str, std::string from, std::string to)
{
	size_t pos = str->find(from, pos);

	while (pos != std::string::npos)
	{
		str->replace(pos, from.length(), to);
		pos += to.length();
		pos = str->find(from, pos);
	}
}

void Parser::outputContents(std::ofstream &stream, int in)
{
	stream << std::setprecision(5);
	stream << indent(in) << "object " << _className
	<< ", " << _absolutePath << std::endl;
	stream << indent(in) << "{" << std::endl;
	in++;

	for (int i = 0; i < _stringProperties.size(); i++)
	{
		std::string name = _stringProperties[i].ptrName;
		std::string *ptr = _stringProperties[i].stringPtr;
		
		sanitise(ptr, "\n", "\\n");
		
		if (!ptr) continue;
		if (ptr->length())        
		{
			stream << indent(in) << name << " = " << *ptr << "" << std::endl;
		}
		else
		{
			stream << indent(in) << name << " = __NULL__" << std::endl;
		}
	}

	for (int i = 0; i < _doubleProperties.size(); i++)
	{
		std::string name = _doubleProperties[i].ptrName;
		double *ptr = _doubleProperties[i].doublePtr;
		if (!ptr) continue;
		stream << indent(in) << name << " = " << f_to_str(*ptr, 10) << std::endl;
	}

	for (int i = 0; i < _mat3x3Properties.size(); i++)
	{
		std::string name = _mat3x3Properties[i].ptrName;
		mat3x3 *ptr = _mat3x3Properties[i].mat3x3Ptr;
		if (!ptr) continue;

		stream << indent(in) << name << " = ";

		for (int i = 0; i < 8; i++)
		{
			stream << ptr->vals[i] << ",";
		}

		stream << ptr->vals[8] << std::endl;
	}

	for (int i = 0; i < _vec3Properties.size(); i++)
	{
		std::string name = _vec3Properties[i].ptrName;
		vec3 *ptr = _vec3Properties[i].vec3Ptr;
		if (!ptr) continue;
		stream << indent(in) << name << " = " << ptr->x
		<< "," << ptr->y << "," << ptr->z << std::endl;
	}

	for (int i = 0; i < _boolProperties.size(); i++)
	{
		std::string name = _boolProperties[i].ptrName;
		bool *ptr = _boolProperties[i].boolPtr;
		if (!ptr) continue;
		stream << indent(in) << name << " = " << *ptr << std::endl;
	}

	for (int i = 0; i < _intProperties.size(); i++)
	{
		std::string name = _intProperties[i].ptrName;
		int *ptr = _intProperties[i].intPtr;
		if (!ptr) continue;
		stream << indent(in) << name << " = " << *ptr << std::endl;
	}

	for (ParserList::iterator it = _parserList.begin();
	     it != _parserList.end(); it++)
	{
		std::string category = it->first;
		stream << indent(in) << "category " << category << std::endl;
		stream << indent(in) << "{" << std::endl;
		in++;       

		for (int i = 0; i < it->second.size(); i++)
		{
			ParserPtr child = it->second.at(i);
			if (!child)
			{
				continue;
			}

			child->outputContents(stream, in);            
		}

		in--;
		stream << indent(in) << "}" << std::endl;
	}

	for (ReferenceList::iterator it = _referenceList.begin();
	     it != _referenceList.end(); it++)
	{
		std::string category = it->first;
		stream << indent(in) << "references " << category << std::endl;
		stream << indent(in) << "{" << std::endl;
		in++;

		for (int j = 0; j < it->second.size(); j++)
		{
			ParserPtr child = it->second.at(j).lock();
			stream << indent(in) << child->getAbsolutePath() << std::endl;
		}

		in--;
		stream << indent(in) << "}" << std::endl;
	}

	for (int i = 0; i < _vec3ArrayProperties.size(); i++)
	{
		Vec3ArrayProperty property = _vec3ArrayProperties[i];
		stream << indent(in) << "array " << property.ptrName << std::endl;
		stream << indent(in) << "{" << std::endl;
		in++;

		for (int j = 0; j < property.vec3ArrayPtr->size(); j++)
		{
			vec3 vec = property.vec3ArrayPtr->at(j);
			stream << indent(in) << vec3_desc(vec) << std::endl;
		}

		in--;
		stream << indent(in) << "}" << std::endl;
	}

	for (int i = 0; i < _mat3x3ArrayProperties.size(); i++)
	{
		Mat3x3ArrayProperty property = _mat3x3ArrayProperties[i];
		stream << indent(in) << "array " << property.ptrName << std::endl;
		stream << indent(in) << "{" << std::endl;
		in++;

		for (int j = 0; j < property.mat3x3ArrayPtr->size(); j++)
		{
			mat3x3 mat = property.mat3x3ArrayPtr->at(j);
			stream << indent(in) << "("; 

			for (int k = 0; k < 8; k++)
			{
				stream << mat.vals[k] << ", ";                
			}

			stream << mat.vals[8] << ")" << std::endl;
		}

		in--;
		stream << indent(in) << "}" << std::endl;
	}

	in--;
	stream << indent(in) << "}" << std::endl;
}

void Parser::restoreState(int num)
{
	int count = _allParsers.size();
	
	if (count == 0) return;
	
	int totalStates = stateCount();
	
	/* If we specify state -1, we want the last-but-one state */
	if (num < 0)
	{
		num = totalStates + num - 1;
	}

	if (num < 0 || num >= totalStates)
	{
		std::cout << "Cannot restore state: " << num << std::endl;
		std::cout << "Total states: " << totalStates << std::endl;
		return;
	}

	for (ParserMap::iterator it = _allParsers.begin();
	     it != _allParsers.end(); it++) 
	{
		if (it->second.expired())
		{
			continue;
		}

		it->second.lock()->privateRestoreState(num);	
	}
	
	for (ParserMap::iterator it = _allParsers.begin();
	     it != _allParsers.end(); it++) 
	{
		if (it->second.expired())
		{
			continue;
		}

		it->second.lock()->_restored = false;
	}

	postRestoreState();	
}

void Parser::postRestoreState()
{

}

void Parser::privateRestoreState(int num)
{
	if (stateCount() < num || _restored)
	{
		return;
	}
	
	StateValueList *list = &_states[num];
	
	for (int i = 0; i < list->size(); i++)
	{
		list->at(i).applyToParser(this);
	}
	
	_restored = true;
	
	/* Remove all the states after the one just restored */
	for (int i = num; i < _states.size(); i++)
	{
		_states.erase(_states.begin() + i);
	}
}

void Parser::saveState()
{
	int count = stateCount();
	int aim = count + 1;
	
	for (ParserMap::iterator it = _allParsers.begin();
	     it != _allParsers.end(); it++) 
	{
		if (it->second.expired())
		{
			continue;
		}

		it->second.lock()->privateSaveState(aim);	
		count = it->second.lock()->stateCount();
	}

	if (count > 100)
	{
		for (ParserMap::iterator it = _allParsers.begin();
		     it != _allParsers.end(); it++) 
		{
			if (it->second.expired())
			{
				continue;
			}

			ParserPtr parser = it->second.lock();
			parser->_states.erase(parser->_states.begin());
		}

	}
}

void Parser::privateSaveState(int aim)
{
	if (stateCount() >= aim)
	{
		return;
	}

	StateValueList list;

	for (int i = 0; i < _stringProperties.size(); i++)
	{
		StateValue value;
		value.addStringValue(_stringProperties[i].ptrName,
		                     *_stringProperties[i].stringPtr);	
		list.push_back(value);
	}

	for (int i = 0; i < _doubleProperties.size(); i++)
	{
		StateValue value;
		value.addDoubleValue(_doubleProperties[i].ptrName,
		                     *_doubleProperties[i].doublePtr);	
		list.push_back(value);
	}

	for (int i = 0; i < _mat3x3Properties.size(); i++)
	{
		StateValue value;
		value.addMat3x3Value(_mat3x3Properties[i].ptrName,
		                     *_mat3x3Properties[i].mat3x3Ptr);	
		list.push_back(value);
	}

	for (int i = 0; i < _vec3Properties.size(); i++)
	{
		StateValue value;
		value.addVec3Value(_vec3Properties[i].ptrName,
		                   *_vec3Properties[i].vec3Ptr);	
		list.push_back(value);
	}

	for (int i = 0; i < _vec3ArrayProperties.size(); i++)
	{
		StateValue value;
		value.addVec3ArrayValue(_vec3ArrayProperties[i].ptrName,
		                        *_vec3ArrayProperties[i].vec3ArrayPtr);	
		list.push_back(value);
	}

	for (int i = 0; i < _mat3x3ArrayProperties.size(); i++)
	{
		StateValue value;
		value.addMat3x3ArrayValue(_mat3x3ArrayProperties[i].ptrName,
		                          *_mat3x3ArrayProperties[i].mat3x3ArrayPtr);
		list.push_back(value);
	}

	for (int i = 0; i < _boolProperties.size(); i++)
	{
		StateValue value;
		value.addBoolValue(_boolProperties[i].ptrName,
		                       *_boolProperties[i].boolPtr);	
		list.push_back(value);
	}

	for (int i = 0; i < _intProperties.size(); i++)
	{
		StateValue value;
		value.addIntValue(_intProperties[i].ptrName,
		                      *_intProperties[i].intPtr);	
		list.push_back(value);
	}

	_states.push_back(list);
}

void Parser::clearContents()
{
	_setup = false;
	_stringProperties.clear();
	_doubleProperties.clear();
	_intProperties.clear();
	_vec3Properties.clear();
	_boolProperties.clear();
	_referenceList.clear();
	_vec3ArrayProperties.clear();
	_mat3x3ArrayProperties.clear();
	_mat3x3Properties.clear();

	for (ParserList::iterator it = _parserList.begin();
	     it != _parserList.end(); it++)
	{
		for (int i = 0; i < it->second.size(); i++)
		{
			ParserPtr child = it->second.at(i);
			if (!child) continue;

			child->clearContents();
		}
	}

	_parserList.clear();
}

void Parser::writeToFile(std::ofstream &stream, int in)
{
	clearContents();
	
	setup();

	stream << "vagabond data structure v0.0" << std::endl;

	outputContents(stream, in);

	stream << std::endl << "END" << std::endl;
}


void Parser::addReference(std::string category, ParserPtr cousin)
{
	if (!cousin) return;
	_referenceList[category].push_back(cousin);
}

char *Parser::parseNextSpecial(char *block)
{
	// posied at the word just after 'special'.

	char *white = strchrwhite(block);
	*white = 0;

	std::string specialName = std::string(block);
	block = white + 1;
	incrementIndent(&block);

	if (block[0] != '{')
	{
		std::cout << "Something's wrong - was expecting a {!" << std::endl;
		return white;
	}

	block++;
	incrementIndent(&block);

	if (block[0] != '}')
	{
		std::cout << "Why has this special thing not ended?" << std::endl;
		std::cout << "Actual char: " << block[0] << std::endl;
	}

	block++;
	incrementIndent(&block);

	return block;
}

char *getCategory(char *block, std::string *catName)
{
	char *white = strchrwhite(block);
	*white = 0;

	std::string categoryName = std::string(block);
//	std::cout << "Category name is " << categoryName << std::endl;
	block = white + 1;
	incrementIndent(&block);

	// now we want an open bracket.

	if (block[0] != '{')
	{
		std::cout << "Something's wrong - was expecting a {!" << std::endl;
		return white;
	}

	block++;
	incrementIndent(&block);

	*catName = categoryName;

	return block;
}

char *Parser::parseNextArray(char *block)
{
	std::string categoryName;
	// poised at the word just after 'array'.
	block = getCategory(block, &categoryName);

	// now poised at whatever came after the {.
	// Now we expect a bundle of ()s.

	// we will start by loading them in as a contiguous array.
	// we shall cast them to the appropriate variable later.
	std::vector<double> stream;

	while (true)
	{
		if (block[0] == '}')
		{
			block++;
			incrementIndent(&block);
			break;
		}

		std::vector<double> values;

		// check that we start with a (.
		if (block[0] == '(')
		{
			block++;
			incrementIndent(&block);
		} 

		bool more = true;

		// now we have a load of numbers separated by commas.
		while (more)
		{
			char *comma = strchr(block, ',');
			char *bracket = strchr(block, ')');

			if (bracket == NULL && comma == NULL)
			{
				std::cout << "Error - truncated file?" << std::endl;
				return NULL;	
			}

			if ((bracket != NULL && bracket < comma) || comma == NULL)
			{
				comma = bracket;
				more = false;
			}

			*comma = 0;

			// we add this to our contiguous stream.  

			double anotherNum = strtod(block, NULL);
			stream.push_back(anotherNum);

			// bring block to the next 'thing'.
			block = comma + 1;
			incrementIndent(&block);
		}

		// eventually we hit a } and we break.
	}

	// find the appropriate array
	for (int i = 0; i < _mat3x3ArrayProperties.size(); i++)
	{
		Mat3x3ArrayProperty property = _mat3x3ArrayProperties[i];
		if (property.ptrName == categoryName)
		{
			size_t size = stream.size() / 9;
			property.mat3x3ArrayPtr->resize(size);
			mat3x3 *start = &((*property.mat3x3ArrayPtr)[0]); 
			memcpy(start, &stream[0], stream.size() * sizeof(double));
		}
	}

	for (int i = 0; i < _vec3ArrayProperties.size(); i++)
	{
		Vec3ArrayProperty property = _vec3ArrayProperties[i];
		if (property.ptrName == categoryName)
		{
			size_t size = stream.size() / 3;
			property.vec3ArrayPtr->resize(size);
			vec3 *start = &((*property.vec3ArrayPtr)[0]);
			memcpy(start, &stream[0], stream.size() * sizeof(double));
		}
	}

	return block;
}

char *Parser::parseNextReference(char *block)
{
	std::string categoryName;
	// poised at the word just after 'references'.
	block = getCategory(block, &categoryName);

	// expecting a fuckton of references now
	// will go through later and fill them in

	while (true)
	{
		char *white = strchrwhite(block);
		*white = 0;

		char *reference = block;

		if (strlen(reference) > 0)
		{
			std::string refStr = std::string(reference);
			_resolveList[categoryName].push_back(refStr);
			//            std::cout << "Adding reference " << reference << std::endl;
		}

		block = white + 1; 
		incrementIndent(&block);

		if (block[0] == '}')
		{
			//            std::cout << "Found a } in references" << std::endl;
			block++;
			incrementIndent(&block);
			return block;
		}

		if (block[0] == 0)
		{
			return NULL;
		}
	}
}

char *Parser::parseNextObject(char *block)
{
	// poised at the word just after 'category'.
	char *white = strchrwhite(block);
	*white = 0;

	std::string categoryName = std::string(block);
	//    std::cout << "Category name is " << categoryName << std::endl;
	block = white + 1;
	incrementIndent(&block);

	// now we want an open bracket.

	if (block[0] != '{')
	{
		std::cout << "Something's wrong - was expecting a {!" << std::endl;
		return white;
	}

	block++;
	incrementIndent(&block);

	// expecting a fuckton of objects now

	white = strchrwhite(block);
	*white = 0;

	if (strncmp(block, "object", 3) != 0)
	{
		std::cout << "Something's wrong - was expecting an object!" << std::endl;
		return white;
	}

	block = white + 1;
	incrementIndent(&block);

	bool stillObjects = true;

	while (stillObjects)
	{
		char *comma = strchr(block, ',');
		*comma = 0;

		ParserPtr object = objectOfType(block);
		if (!object)
		{
			return NULL;
		}

		block = comma + 1;

		// Comes back incremented.
		block = object->parse(block);

		addObject(object, categoryName);
		_parserList[categoryName].push_back(object);
		std::string path = object->getAbsolutePath();

		if (_allParsers.count(path) > 0)
		{
			shout_at_user("Duplicate object " + path);
		}

		addToAllParsers(path, object);

		if (block == NULL)
		{
			return NULL;
		}

		if (block[0] == 0)
		{
			std::cout << "File ended prematurely!" << std::endl;
			return NULL;
		}

		white = strchrwhite(block);
		*white = 0;

		if (strncmp(block, "object", 3) != 0)
		{
			stillObjects = false;
		}
		else
		{
			block = white + 1;
		}
	}

	if (block[0] == '}')
	{
		block++;
		incrementIndent(&block);
	}

	return block;
}

void Parser::setProperty(std::string property, std::string value)
{
	for (int i = 0; i < _stringProperties.size(); i++)
	{
		if (_stringProperties[i].ptrName == property)
		{
			if (value == "__NULL__")
			{
				*_stringProperties[i].stringPtr = "";
				return; 
			}

			sanitise(&value, "\\n", "\n");
			*_stringProperties[i].stringPtr = value;

			return;
		}
	}

	for (int i = 0; i < _doubleProperties.size(); i++)
	{
		if (_doubleProperties[i].ptrName == property)
		{
			char *check = NULL;
			double val = strtod(value.c_str(), &check);            
			if (check > &value[0]);
			{
				*_doubleProperties[i].doublePtr = val;
			}
			return;
		}
	}

	for (int i = 0; i < _mat3x3Properties.size(); i++)
	{
		if (_mat3x3Properties[i].ptrName == property)
		{
			mat3x3 newMat = make_mat3x3();
			char *start = &value[0];

			for (int i = 0; i < 9; i++)
			{
				char *comma = strchr(start, ',');

				// last one needs no treatment, null-terminated already
				if (i < 8)
				{
					*comma = 0;
				}

				char *check = NULL;
				double x = strtod(start, &check);

				if (check == NULL)
				{
					std::cout << "Error while parsing matrix property in "
					<< getAbsolutePath() << std::endl;
					break;
				}

				newMat.vals[i] = x;
				start = comma + 1;
			}

			*_mat3x3Properties[i].mat3x3Ptr = newMat;
			return;
		}
	}

	for (int i = 0; i < _vec3Properties.size(); i++)
	{
		if (_vec3Properties[i].ptrName == property)
		{
			char *start = &value[0];
			char *comma = strchr(start, ',');
			*comma = 0;
			double x = strtod(start, NULL);
			start = comma + 1;
			comma = strchr(start, ',');
			*comma = 0;
			double y = strtod(start, NULL);
			start = comma + 1;

			// Next one is already zero 
			double z = strtod(start, NULL);
			vec3 vec = make_vec3(x, y, z); 
			*_vec3Properties[i].vec3Ptr = vec;
			return;
		}
	}

	for (int i = 0; i < _intProperties.size(); i++)
	{
		if (_intProperties[i].ptrName == property)
		{
			int val = atoi(value.c_str());
			*_intProperties[i].intPtr = val;
			return;
		}
	}

	for (int i = 0; i < _boolProperties.size(); i++)
	{
		if (_boolProperties[i].ptrName == property)
		{
			bool val = atoi(value.c_str());
			*_boolProperties[i].boolPtr = val;
			return;
		}
	}

	std::cout << "Unhandled thing: " << property << " = " << value << "." << std::endl;
}

char *Parser::parseNextProperty(std::string property, char *block)
{
	// just incremented to = sign, let's check.
	if (block[0] != '=')
	{
		std::cout << "Failure to read property for " << _identifier << ", no = sign! It's " << block[0] << std::endl;
		return &block[1];
	}

	block++;
	incrementIndent(&block);
	char *white = strchr(block, '\n');
	*white = '\0';

	// now we expect the value between block and white...
	std::string value = std::string(block);

	setProperty(property, value);

	block = white + 1;
	incrementIndent(&block);

	return block;
}

bool Parser::parseNextChunk(char **blockPtr)
{
	char *block = *blockPtr;    

	if (block == NULL)
	{
		std::cout << "On object " << getAbsolutePath() << std::endl;
		std::cout << "Parsing error has occurred." << std::endl;
		return false;
	}

	// just incremented from the first {.
	// we should expect a keyword, or a property next.
	char *space = strchrwhite(block);

	if (space == NULL)
	{
		std::cout << "Nope!" << std::endl;
		std::cout << block << std::endl;
		return false;
	}

	*space = '\0';

	// what are we about to process?

	ParserType type = ParserTypeProperty;

	if (strncmp(block, "category", 4) == 0)
	{
		type = ParserTypeObject;
	}
	else if (strncmp(block, "references", 4) == 0)
	{
		type = ParserTypeReference;
	}
	else if (strncmp(block, "array", 4) == 0)
	{
		type = ParserTypeArray;
	}
	else if (strncmp(block, "special", 4) == 0)
	{
		type = ParserTypeSpecial;
	}

	char *property = block;
	block = space + 1;
	incrementIndent(&block);

	switch (type)
	{
		case ParserTypeProperty:
		*blockPtr = parseNextProperty(property, block);
		break;

		case ParserTypeObject:
		*blockPtr = parseNextObject(block);
		break;

		case ParserTypeSpecial:
		*blockPtr = parseNextSpecial(block);
		break;

		case ParserTypeReference:
		*blockPtr = parseNextReference(block);
		break;

		case ParserTypeArray:
		*blockPtr = parseNextArray(block);
		break;

		default:
		break;
	}

	block = *blockPtr;

	if (*blockPtr == NULL)
	{
		std::cout << "Parsing error occurred." << std::endl;
		return false;
	}

	if (block[0] == '}')
	{
		block++;
		incrementIndent(&block);
		return false;
	}

	return true;
}

char *Parser::parse(char *block)
{
	// Get to the beginning of the absolute path...
	incrementIndent(&block);

	char *newline = strchrwhite(block);
	if (newline == NULL) return NULL;

	// We prepare the vessels to accept our data.
	setup(true);

	// Now we can replace the dud path/identifier.
	*newline = '\0';
	std::string path = std::string(block);
	_absolutePath = path;

	// we also want to isolate the final identifier
	char *slash = strrchr(block, '/');
	slash++;
	std::string identifier = std::string(slash);
	_identifier = identifier;

	//    std::cout << "Found identifier " << _identifier << std::endl;

	block = newline + 1;
	incrementIndent(&block);
	// we expect the next character to be a {
	if (block[0] != '{')
	{
		std::cout << "Check me: no { for object " << path  << "?" << std::endl;
		return NULL;
	}

	// Get past the {. 
	block++;
	incrementIndent(&block);

	// Loop through the properties.

	bool another = true;

	while (another)
	{    
		another = parseNextChunk(&block);
	}

	if (block == NULL)
	{
		return NULL;
	}

	if (block[0] == 0)
	{
		return NULL;
	}

	if (block[0] != '}')
	{
		std::cout << "Why is there no }?" << std::endl;    
	}
	//    else std::cout << "There is an } at the end of an object" << std::endl;    

	block++;
	incrementIndent(&block);

	return block;
}

ParserPtr Parser::objectOfType(char *className)
{
	ParserPtr object = ParserPtr();

	if (strcmp(className, "Crystal") == 0) 
	{
		object = ParserPtr(static_cast<Crystal *>(new Crystal()));
		Options::getRuntimeOptions()->addCrystal(object);
	}
	else if (strcmp(className, "Polymer") == 0)
	{
		object = ParserPtr(static_cast<Polymer *>(new Polymer()));
	}
	else if (strcmp(className, "Atom") == 0)
	{
		object = ParserPtr(static_cast<Atom *>(new Atom()));
	}
	else if (strcmp(className, "Absolute") == 0)
	{
		object = ParserPtr(static_cast<Absolute *>(new Absolute()));        
	}
	else if (strcmp(className, "Anchor") == 0)
	{
		object = ParserPtr(static_cast<Anchor *>(new Anchor()));        
	}
	else if (strcmp(className, "Bond") == 0)
	{
		object = ParserPtr(static_cast<Bond *>(new Bond()));        
	}
	else if (strcmp(className, "Ghost") == 0)
	{
		object = ParserPtr(static_cast<GhostBond *>(new GhostBond()));        
	}
	else if (strcmp(className, "Whack") == 0)
	{
		object = ParserPtr(static_cast<Whack *>(new Whack()));        
	}
	else if (strcmp(className, "BondGroup") == 0)
	{
		object = ParserPtr(static_cast<BondGroup *>(new BondGroup()));        
	}
	else if (strcmp(className, "Molecule") == 0)
	{
		object = ParserPtr(static_cast<Molecule *>(new Molecule()));        
	}
	else if (strcmp(className, "WaterNetwork") == 0)
	{
		object = ParserPtr(static_cast<WaterNetwork *>(new WaterNetwork()));        
	}
	else if (strcmp(className, "Monomer") == 0)
	{
		object = ParserPtr(static_cast<Monomer *>(new Monomer()));        
	}
	else if (strcmp(className, "Sidechain") == 0)
	{
		object = ParserPtr(static_cast<Sidechain *>(new Sidechain()));        
	}
	else if (strcmp(className, "Backbone") == 0)
	{
		object = ParserPtr(static_cast<Backbone *>(new Backbone()));        
	}
	else
	{
		std::cout << "Do not understand class name " << className << std::endl;
		return object;
	}

	//    std::cout << "Making object of type " << className << std::endl;
	return object;
}

void Parser::setupKnownClasses()
{
	_allClasses["Molecule"] = 1;
	_allClasses["Model"] = 1;
}

void Parser::addToAllParsers(std::string key, ParserPtr parser)
{
	_allParsers[key] = parser;
	std::string name = parser->getClassName();

	if (_allClasses.count(name) == 0)
	{
		_allClasses[name] = 0;
	}
	
	_allClasses[name]++;
}

ParserPtr Parser::processBlock(char *block)
{
	_allParsers.clear();

	char *comma = strchr(block, ',');
	*comma = '\0';

	// Get the initial object (crystal, we hope) 
	ParserPtr object = objectOfType(block);
	if (!object)
	{
		return ParserPtr();
	}

	block = comma + 1;

	// Parse the entirety of the structure.
	char *success = object->parse(block);

	std::cout << "Vagabond file parsing complete..." << std::endl;
	// Resolve dangling references.
	object->resolveReferences();

	std::cout << "Object references resolved..." << std::endl;
	// Add parent to complete parser list
	addToAllParsers(object->getAbsolutePath(), object);

	// Loop through all objects to allow them to finish up.
	for (ParserMap::iterator it = _allParsers.begin();
	     it != _allParsers.end(); it++)
	{
		ParserPtr aParser = it->second.lock();
		aParser->postParseTidy();
	}


	object->postParseTidy();

	std::cout << "Post-parse object tidy done..." << std::endl;

	if (success != NULL)
	{
		return object;
	}

	return ParserPtr();
}

ParserPtr Parser::resolveReference(std::string reference)
{
	ParserMap::iterator it;
	it = _allParsers.find(reference);
	
	if (it == _allParsers.end())
	{
		return ParserPtr();
	}
	
	return it->second.lock();
}

void Parser::resolveReferences()
{
	for (ResolveList::iterator it = _resolveList.begin();
	     it != _resolveList.end(); it++)
	{
		for (int j = 0; j < it->second.size(); j++)
		{
			std::string path = it->second[j];
			ParserPtr parser = resolveReference(path);
			linkReference(parser, it->first);
		}
	}

	for (ParserList::iterator it = _parserList.begin();
	     it != _parserList.end(); it++) 
	{
		for (int j = 0; j < it->second.size(); j++)
		{
			ParserPtr child = it->second[j];
			if (!child) continue;

			child->resolveReferences();
		}
	}
}

