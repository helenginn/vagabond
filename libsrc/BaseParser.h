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
#include <hcsrc/RefinementStrategy.h>

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
 * linkReference(BaseParserPtr, std::string) needs to be added. In some cases,
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
typedef std::map<std::string, TwoDouble > TwoDoubleList;
typedef std::map<std::string, BaseParserWkr> ParserMap;
typedef std::map<std::string, int> ClassMap;

class StateValue;
typedef std::vector<StateValue> StateValueList;

class BaseParser : public boost::enable_shared_from_this<BaseParser>
{
public:
	friend class StateValue;
	BaseParser();
	virtual ~BaseParser() {};
	
	void prepare();

	std::string getAbsolutePath()
	{
		return _absolutePath;
	}

	/**
	* 	Implementation of getClassName in downstream classes should 
	*   return the name of the class itself.
	*/
	virtual std::string getClassName() = 0;
	
	/**
	 * Outside of the class this should only be called on the top level of
	 * the contents of a Vagabond file, and will return a BaseParserPtr of subclass
	 * Crystal. */
	static BaseParserPtr processBlock(char *block);

	/** \return returns number of saved 'undo' states available */
	size_t stateCount()
	{
		return _states.size();	
	}
	
	/** Save a state which can then be subject to 'undo' */
	void saveState();
	
	/** Restore a previous state.
	 * \param num if positive, will take the Nth saved state; if negative,
	 * will take the -Nth from the end. */
	void restoreState(int num);
	
	virtual std::string shortDesc()
	{
		return "";
	}
	
	virtual std::string longDesc()
	{
		return "";
	}

	/** Number of classes in the current structure, to be called on a top
	 * level Crystal BaseParserPtr. */
	int classCount(std::string name)
	{
		if (_top->_allClasses.count(name) == 0)
		{
			return 0;
		}
		
		return _allClasses[name];
	}
	
	/**
	*  Once all the cross-references to Parsers have been completed and all
	*  Parsers exist in memory, this function is called. Here it may be
	*  appropriate to perform some final tidying functions.	
	*/
	virtual void postParseTidy() {};
	void postParse();
	
	void setTop(BaseParserPtr top)
	{
		_top = &*top;
	}
	
	BaseParserPtr getTopParser()
	{
		return _top->shared_from_this();
	}
	
	void outputStream();

	virtual void setStream(std::ostream *str)
	{
		_stream = str;
	}

	std::ostream *getStream()
	{
		return _stream;
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
	* 	and addMat3x3Property(). Also use addChild() here
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
	virtual void linkReference(BaseParserPtr, std::string) {};
	
	/** \name Adding property functions */
	/**@{*/
	
	/**
	 * Add a reference to the pointer to a std::string. Should not be one
	 * which supplies line breaks as this will mess up the file format - no
	 * escaping done at the moment. */
	void addStringProperty(std::string className, std::string *ptr);
	
	/**
	 * Adds a reference to a pointer to a double. Will be reproduced to
	 * 10 decimal places. */
	void addDoubleProperty(std::string className, double *ptr);

	/**
	 * Adds a reference to a pointer to an integer. */
	void addIntProperty(std::string className, int *ptr);
	
	/**
	 * Adds a reference to a vec3 pointer */
	void addVec3Property(std::string className, vec3 *ptr);

	/**
	 * Adds a reference to a mat3x3 pointer */
	void addMat3x3Property(std::string className, mat3x3 *ptr, 
	                       bool onlyReceive = false);
	void addBoolProperty(std::string className, bool *ptr);
	void addChild(std::string category, ParserPtr child);
	void addReference(std::string category, ParserPtr cousin);
	void addVec3ArrayProperty(std::string className, std::vector<vec3> *ptr,
	                          bool receiveOnly = false);
	void addIntArrayProperty(std::string className, 
	                         std::vector<unsigned long> *ptr,
	                          bool receiveOnly = false);
	void addMat3x3ArrayProperty(std::string className, std::vector<mat3x3> *ptr);
	/**@}*/


	bool hasSetter(std::string funcName)
	{
		return (_setterList.count(funcName) > 0);		
	}
	
	bool hasTwoDouble(std::string funcName)
	{
		return (_twoDoubleList.count(funcName) > 0);		
	}
	
	TwoDouble getTwoDouble(std::string funcName)
	{
		return _twoDoubleList[funcName];
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

	BaseParserPtr resolveReference(std::string reference);
	
	void privateRestoreState(int num);
	virtual void postRestoreState();

	std::ostream *_stream;
private:
	std::string _className;
	std::string _identifier;
	std::string _absolutePath;    
	BaseParser *_parent;
	BaseParser *_top;
	std::vector<StateValueList> _states;

	std::vector<StringProperty> _stringProperties;
	std::vector<DoubleProperty> _doubleProperties;
	std::vector<IntProperty> _intProperties;
	std::vector<Vec3Property> _vec3Properties;
	std::vector<Mat3x3Property> _mat3x3Properties;
	std::vector<BoolProperty> _boolProperties;
	std::vector<Vec3ArrayProperty> _vec3ArrayProperties;
	std::vector<IntArrayProperty> _intArrayProperties;
	std::vector<Mat3x3ArrayProperty> _mat3x3ArrayProperties;
	ResolveList _resolveList;
	ParserList _parserList;
	ReferenceList _referenceList;
	FunctionList _functionList;
	TwoDoubleList _twoDoubleList;
	SetterList _setterList;

	void makePath();
	void setup(bool isNew = false);
	bool _setup;
	bool _restored;
	
	void privateSaveState(int aim);
	
	BaseParser *getParent()
	{
		return _parent;	
	}

	void addToAllParsers(std::string key, BaseParserPtr parser);
	void outputContents(std::ofstream &stream, int in);
	void clearContents();
	void setParent(BaseParser *parent);
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

	ParserMap _allParsers;
	ClassMap _allClasses;
};


#endif
