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
#include <sstream>
#include <iostream>

typedef enum
{
    ParserTypeObject,
    ParserTypeReference,
    ParserTypeSpecial,
    ParserTypeProperty,
} ParserType;

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
    bool *boolPtr; 
    std::string ptrName;
} BoolProperty;

typedef struct
{
    int *intPtr;
    std::string ptrName;
} IntProperty;

typedef void (*Encoder)(void *, void *, std::ofstream &, int);
typedef char *(*Decoder)(void *, void *, char *block);

typedef struct
{
    void *objPtr;
    std::string ptrName;
    void *delegate;
    Encoder encoder;
    Decoder decoder;
} CustomProperty;

inline std::string indent(int num)
{
    std::ostringstream stream;

    for (int i = 0; i < num; i++)
    {
        stream << "  ";
    }
    return stream.str();
}

inline char *strchrwhite(char *block)
{
    char *space = strchr(block, ' ');
    char *newline = strchr(block, '\n');
    char *tab = strchr(block, '\t');

    if (tab != NULL && tab < newline) newline = tab;
    if (space != NULL && space < newline) newline = space;

    return newline;
}

inline void incrementIndent(char **block)
{
    while ((*block)[0] == ' ' || (*block)[0] == '\t' || (*block)[0] == '\n'
            || (*block)[0] == '\0')
    {
        (*block)++;
    }
}

inline char *keywordValue(char *block, char **keyword, char **value) 
{
    char *space = strchrwhite(block);

    if (space == NULL)
    {
        std::cout << "Space is just null" << std::endl;
        return NULL;
    }

    *space = '\0';
    *keyword = block;
    block = space + 1;
    incrementIndent(&block);
    
    // Don't panic, we probably just have an 'object'.
    if (block[0] != '=')
    {
        std::cout << "keyword: " << *keyword << " - block char " << *block << std::endl;
        return block;
    }

    block++;
    incrementIndent(&block);
    
    space = strchrwhite(block);
    *space = '\0';
    
    *value = block;
    block = space + 1;
    incrementIndent(&block);

    return block;
}


typedef std::map<std::string, std::vector<ParserPtr> > ParserList;
typedef std::map<std::string, std::vector<ParserPtr> > ReferenceList;
typedef std::map<std::string, std::vector<std::string> > ResolveList;

class Parser 
{
public:
    Parser();

    std::string getAbsolutePath()
    {
        return _absolutePath;
    }

    virtual std::string getClassName() = 0;
    static ParserPtr processBlock(char *block);
protected:
    virtual std::string getParserIdentifier() = 0;
    virtual void addProperties() = 0;
    virtual void addObject(ParserPtr object, std::string category) {};

    void setParent(Parser *parent);
    void addStringProperty(std::string className, std::string *ptr);
    void addDoubleProperty(std::string className, double *ptr);
    void addIntProperty(std::string className, int *ptr);
    void addVec3Property(std::string className, vec3 *ptr);
    void addCustomProperty(std::string className, void *ptr,
                           void *delegate, Encoder encoder,
                           Decoder decoder); 
    void addBoolProperty(std::string className, bool *ptr);
    void addChild(std::string category, ParserPtr child);
    void addReference(std::string category, ParserPtr cousin);

    void writeToFile(std::ofstream &stream, int indent);
    void clearContents();

private:
    std::string _className;
    std::string _identifier;
    std::string _absolutePath;    
    Parser *_parent;

    std::vector<StringProperty> _stringProperties;
    std::vector<DoubleProperty> _doubleProperties;
    std::vector<IntProperty> _intProperties;
    std::vector<Vec3Property> _vec3Properties;
    std::vector<BoolProperty> _boolProperties;
    std::vector<CustomProperty> _customProperties;
    ResolveList _resolveList;
    ParserList _parserList;
    ReferenceList _referenceList;

    void makePath();
    void setup(bool isNew = false);
    bool _setup;
    
    void outputContents(std::ofstream &stream, int in);
    static ParserPtr objectOfType(char *className);
    bool parseNextChunk(char **blockPtr);
    char *parse(char *block);
    char *parseNextProperty(std::string property, char *block);
    char *parseNextObject(char *block);
    char *parseNextSpecial(char *block);
    char *parseNextReference(char *block);
    void setProperty(std::string property, std::string value);

    static ParserList _allParsers;
};


#endif
