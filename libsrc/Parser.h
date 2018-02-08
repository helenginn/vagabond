// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#ifndef __vagabond__Parser__
#define __vagabond__Parser__

#include "shared_ptrs.h"
#include <map>
#include <fstream>
#include "vec3.h"

typedef struct
{
    std::string *stringPtr;
    std::string ptrName;
} StringProperty;

typedef struct
{
    double *doublePtr;
    std::string ptrName;
} DoubleProperty;

typedef struct
{
    vec3 *vec3Ptr;
    std::string ptrName;
} Vec3Property;

typedef struct
{
    int *intPtr;
    std::string ptrName;
} IntProperty;

typedef void (*Encoder)(void *, void *, std::ofstream &, int);

typedef struct
{
    void *objPtr;
    std::string ptrName;
    void *delegate;
    Encoder encoder;
} CustomProperty;

typedef std::map<std::string, std::vector<ParserPtr> > ParserList;
typedef std::map<std::string, std::vector<ParserPtr> > ReferenceList;

class Parser 
{
public:
    Parser();

protected:
    virtual std::string getClassName() = 0;
    virtual std::string getParserIdentifier() = 0;
    virtual void addProperties() = 0;

    void setParent(Parser *parent);
    void addStringProperty(std::string className, std::string *ptr);
    void addDoubleProperty(std::string className, double *ptr);
    void addIntProperty(std::string className, int *ptr);
    void addVec3Property(std::string className, vec3 *ptr);
    void addCustomProperty(std::string className, void *ptr,
                           void *delegate, Encoder encoder); 
    void addChild(std::string category, ParserPtr child);
    void addReference(std::string category, ParserPtr cousin);

    void writeToFile(std::ofstream &stream, int indent);
private:
    std::string _className;
    std::string _identifier;
    std::string _absolutePath;    
    Parser *_parent;

    std::string getAbsolutePath()
    {
        return _absolutePath;
    }

    std::vector<StringProperty> _stringProperties;
    std::vector<DoubleProperty> _doubleProperties;
    std::vector<IntProperty> _intProperties;
    std::vector<Vec3Property> _vec3Properties;
    std::vector<CustomProperty> _customProperties;
    ParserList _parserList;
    ReferenceList _referenceList;

    void makePath();
    void setup();
    bool _setup;
    
    void outputContents(std::ofstream &stream, int in);
};


#endif
