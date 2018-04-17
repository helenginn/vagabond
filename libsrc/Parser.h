// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#ifndef __vagabond__Parser__
#define __vagabond__Parser__

#include "ParserTypes.h"
#include "shared_ptrs.h"
#include <map>
#include <fstream>
#include <iostream>
#include "StateValue.h"
#include "RefinementStrategy.h"

/**
 * \class Parser
 * \brief Parser objects provide support for reading and writing from Vagabond
 * file formats.
 * 
 * An Object of a type Parser can be organised in a hierarchy which is then
 * read and written to disk in a .vbond format. To be a Parser, a small number
 * of functions must be implemented in subclasses.
 *
 * The three functions getClassName(), addProperties() and
 * getParserIdentifier() must be implemented for writing to files. Use
 * addChild() to add other Parsers in a hierarchical fashion.
 *
 * To read back, basic properties (doubles, ints, vectors of int, vec3, mat3x3
 * etc.) are automatically populated. Objects are given back to a Parser from
 * a file using addObject(). If child objects have been specified,
 * linkReference(ParserPtr, std::string) needs to be added. In some cases,
 * postParseTidy() might need to be implemented. This is only called once all
 * Parser objects have been initialised and are in place.
 */

typedef enum
{
	ParserTypeObject,
	ParserTypeReference,
	ParserTypeSpecial,
	ParserTypeProperty,
	ParserTypeArray,
} ParserType;


typedef std::map<std::string, std::vector<ParserPtr> > ParserList;
typedef std::map<std::string, std::vector<ParserWkr> > ReferenceList;
typedef std::map<std::string, std::vector<std::string> > ResolveList;
typedef std::map<std::string, Getter > FunctionList;
typedef std::map<std::string, Setter > SetterList;
typedef std::map<std::string, ParserWkr> ParserMap;
typedef std::map<std::string, int> ClassMap;

class StateValue;
typedef std::vector<StateValue> StateValueList;

class Parser : public boost::enable_shared_from_this<Parser>
{
public:
	friend class StateValue;
	Parser();
	virtual ~Parser() {};

	std::string getAbsolutePath()
	{
		return _absolutePath;
	}

	/**
	* 	Implementation of getClassName should return the name of the class
	* 	itself.
	*/
	virtual std::string getClassName() = 0;
	static ParserPtr processBlock(char *block);

	size_t stateCount()
	{
		return _states.size();	
	}
	
	static void saveState();
	static void restoreState(int num);

	static int classCount(std::string name)
	{
		if (_allClasses.count(name) == 0)
		{
			return 0;
		}
		
		return _allClasses[name];
	}
protected:
	/**
	* 	Implementation of the parser identifier should return a name of the
	* 	Parser object which is unique within a set of siblings. This will be
	* 	incorporated into a long, descriptive name including the hierarchy.
	*/
	virtual std::string getParserIdentifier() = 0;
	
	/**
	* 	All basic properties should be added in this implementation by calling
	* 	the various property functions such as addStringProperty(),
	* 	addDoubleProperty(), addIntProperty(), addVec3Property(),
	* 	addMat3x3Property() and addCustomProperty(). Also use addChild() here
	* 	to add other Parser objects. */
	virtual void addProperties() = 0;
	
	/**
	*   This function should be implemented to add an object back from a file,
	*   which was specified using addChild() in addProperties().
	*/
	virtual void addObject(ParserPtr /* parser */, std::string /*name */) {};
	
	/**
	* 	This function will be called when all Parser objects have been
	* 	initialised, and now the cross-references to Parsers have to be added.
	*/
	virtual void linkReference(ParserPtr, std::string) {};
	
	/**
	*  Once all the cross-references to Parsers have been completed and all
	*  Parsers exist in memory, this function is called. Here it may be
	*  appropriate to perform some final tidying functions.	
	*/
	virtual void postParseTidy() {};

	/** \name Adding property functions */
	/**@{*/
	void addStringProperty(std::string className, std::string *ptr);
	void addDoubleProperty(std::string className, double *ptr);
	void addIntProperty(std::string className, int *ptr);
	void addVec3Property(std::string className, vec3 *ptr);
	void addMat3x3Property(std::string className, mat3x3 *ptr);
	void addCustomProperty(std::string className, void *ptr,
	                       void *delegate, Encoder encoder,
	Decoder decoder); 
	void addBoolProperty(std::string className, bool *ptr);
	void addChild(std::string category, ParserPtr child);
	void addReference(std::string category, ParserPtr cousin);
	void addVec3ArrayProperty(std::string className, std::vector<vec3> *ptr);
	void addMat3x3ArrayProperty(std::string className, std::vector<mat3x3> *ptr);
	void exposeFunction(std::string funcName, Getter func);
	void exposeFunction(std::string funcName, Setter func);
	/**@}*/

	bool hasFunction(std::string funcName)
	{
		return (_functionList.count(funcName) > 0);		
	}
	
	Getter getFunction(std::string funcName)
	{
		return _functionList[funcName];
	}
	
	bool hasSetter(std::string funcName)
	{
		return (_setterList.count(funcName) > 0);		
	}
	
	Setter getSetter(std::string funcName)
	{
		return _setterList[funcName];
	}
	
friend class Thing;
	std::string *getStringProperty(std::string className);
	double *getDoubleProperty(std::string className);
	int *getIntProperty(std::string className);
	ParserPtr getChild(std::string className, int num);
	int getChildCount(std::string className);
	
	/**
	* Top level object should be called to write to a stream (usually a
	* Crystal).
	*/
	void writeToFile(std::ofstream &stream, int indent);

	static ParserPtr resolveReference(std::string reference);
	
	void privateRestoreState(int num);
	virtual void postRestoreState() {};
private:
	std::string _className;
	std::string _identifier;
	std::string _absolutePath;    
	Parser *_parent;
	std::vector<StateValueList> _states;

	std::vector<StringProperty> _stringProperties;
	std::vector<DoubleProperty> _doubleProperties;
	std::vector<IntProperty> _intProperties;
	std::vector<Vec3Property> _vec3Properties;
	std::vector<Mat3x3Property> _mat3x3Properties;
	std::vector<BoolProperty> _boolProperties;
	std::vector<CustomProperty> _customProperties;
	std::vector<Vec3ArrayProperty> _vec3ArrayProperties;
	std::vector<Mat3x3ArrayProperty> _mat3x3ArrayProperties;
	ResolveList _resolveList;
	ParserList _parserList;
	ReferenceList _referenceList;
	FunctionList _functionList;
	SetterList _setterList;

	void makePath();
	void setup(bool isNew = false);
	bool _setup;
	
	void privateSaveState();
	
	Parser *getParent()
	{
		return _parent;	
	}

	static void addToAllParsers(std::string key, ParserPtr parser);
	void outputContents(std::ofstream &stream, int in);
	void clearContents();
	void setParent(Parser *parent);
	static ParserPtr objectOfType(char *className);
	bool parseNextChunk(char **blockPtr);
	char *parse(char *block);
	char *parseNextProperty(std::string property, char *block);
	char *parseNextObject(char *block);
	char *parseNextSpecial(char *block);
	char *parseNextReference(char *block);
	char *parseNextArray(char *block);
	void setProperty(std::string property, std::string value);
	void resolveReferences();

	static ParserMap _allParsers;
	static ClassMap _allClasses;
};


#endif
