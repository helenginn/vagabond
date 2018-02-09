// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#include "Parser.h"
#include <iostream>
#include <iomanip>
#include "Crystal.h"
#include "Polymer.h"

Parser::Parser()
{
    _setup = false;
    _parent = NULL;
}

void Parser::setup()
{
    if (_setup) return;

    _identifier = getParserIdentifier(); 
    _className = getClassName();

    makePath();
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

void Parser::addDoubleProperty(std::string className, double *ptr)
{
    DoubleProperty property;
    property.ptrName = className;
    property.doublePtr = ptr;
    _doubleProperties.push_back(property);
}

void Parser::addIntProperty(std::string className, int *ptr)
{
    IntProperty property;
    property.ptrName = className;
    property.intPtr = ptr;
    _intProperties.push_back(property);
}

void Parser::addVec3Property(std::string className, vec3 *ptr)
{
    Vec3Property property;
    property.ptrName = className;
    property.vec3Ptr = ptr;
    _vec3Properties.push_back(property);
}

void Parser::addBoolProperty(std::string className, bool *ptr)
{
    BoolProperty property;
    property.ptrName = className;
    property.boolPtr = ptr;
    _boolProperties.push_back(property);
}

void Parser::addCustomProperty(std::string className, void *ptr,
                               void *delegate, Encoder encoder)
{
    CustomProperty property;
    property.ptrName = className;
    property.objPtr = ptr;
    property.delegate = delegate;
    property.encoder = encoder;
    _customProperties.push_back(property);
}

void Parser::addChild(std::string category, ParserPtr child)
{
    if (!child) return;

    child->setParent(this);
    _parserList[category].push_back(child);
    child->setup();
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
        if (!ptr) continue;
        stream << indent(in) << name << " = " << *ptr << "" << std::endl;
    }
 
    for (int i = 0; i < _doubleProperties.size(); i++)
    {
        std::string name = _doubleProperties[i].ptrName;
        double *ptr = _doubleProperties[i].doublePtr;
        if (!ptr) continue;
        stream << indent(in) << name << " = " << *ptr << std::endl;
    }

    for (int i = 0; i < _vec3Properties.size(); i++)
    {
        std::string name = _vec3Properties[i].ptrName;
        vec3 *ptr = _vec3Properties[i].vec3Ptr;
        if (!ptr) continue;
        stream << indent(in) << name << " = " << ptr->x
               << ", " << ptr->y << ", " << ptr->z << std::endl;
    }

    for (int i = 0; i < _customProperties.size(); i++)
    {
        std::string name = _customProperties[i].ptrName;
        void *ptr = _customProperties[i].objPtr;
        void *delegate = _customProperties[i].delegate;
        if (!ptr) continue;

        stream << indent(in) << "special " << name << std::endl;  
        stream << indent(in) << "{" << std::endl;
        in++;
        Encoder encoder = _customProperties[i].encoder;
        (*encoder)(delegate, ptr, stream, in);
        in--;
        stream << indent(in) << "}" << std::endl;
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
            ParserPtr child = it->second.at(j);
            stream << indent(in) << child->getAbsolutePath() << std::endl;
        }

        in--;
        stream << indent(in) << "}" << std::endl;
    }

    in--;
    stream << indent(in) << "}" << std::endl;
}

void Parser::clearContents()
{
    _stringProperties.clear();
    _doubleProperties.clear();
    _intProperties.clear();
    _vec3Properties.clear();
    _customProperties.clear();
    _boolProperties.clear();
    _referenceList.clear();

    for (ParserList::iterator it = _parserList.begin();
            it != _parserList.end(); it++)
    {
        for (int i = 0; i < it->second.size(); i++)
        {
            ParserPtr child = it->second.at(i);
            child->clearContents();
        }
    }

    _parserList.clear();
}

void Parser::writeToFile(std::ofstream &stream, int in)
{
    setup();

    stream << "vagabond data structure v0.0" << std::endl;

    outputContents(stream, in);
    clearContents();
}


void Parser::addReference(std::string category, ParserPtr cousin)
{
    if (!cousin) return;
    _referenceList[category].push_back(cousin);
}

char *strchrwhite(char *block)
{
    char *space = strchr(block, ' ');
    char *newline = strchr(block, '\n');
    char *tab = strchr(block, '\t');

    if (tab != NULL && tab < newline) newline = tab;
    if (space != NULL && space < newline) newline = space;

    return newline;
}

void incrementIndent(char **block)
{
    while ((*block)[0] == ' ' || (*block)[0] == '\t' || (*block)[0] == '\n'
            || (*block)[0] == '\0')
    {
        (*block)++;
    }
}

void Parser::setProperty(std::string property, std::string value)
{
    for (int i = 0; i < _stringProperties.size(); i++)
    {
        if (_stringProperties[i].ptrName == property)
        {
            *_stringProperties[i].stringPtr = value;
            std::cout << "Setting string to " << value << std::endl;
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
            std::cout << "Setting double to " << val << std::endl;
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
            start = comma++;
            incrementIndent(&start);
            comma = strchr(start, ',');
            *comma = 0;
            double y = strtod(start, NULL);
            start = comma++;
            incrementIndent(&start);
            comma = strchrwhite(start);
            *comma = 0;
            double z = strtod(start, NULL);
            vec3 vec = make_vec3(x, y, z); 
            *_vec3Properties[i].vec3Ptr = vec;
            std::cout << "Setting vec3 to " << vec3_desc(vec) << std::endl;
            return;
        }
    }

    for (int i = 0; i < _intProperties.size(); i++)
    {
        if (_intProperties[i].ptrName == property)
        {
            int val = atoi(value.c_str());
            std::cout << "Setting int to " << val << std::endl;
            *_intProperties[i].intPtr = val;
            return;
        }
    }

    for (int i = 0; i < _boolProperties.size(); i++)
    {
        if (_boolProperties[i].ptrName == property)
        {
            bool val = atoi(value.c_str());
            std::cout << "Setting bool to " << val << std::endl;
            *_boolProperties[i].boolPtr = val;
            return;
        }
    }

    std::cout << "Unhandled thing." << std::endl;
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
    char *white = strchrwhite(block); 
    *white = '\0';
    
    // now we expect the value between block and white...
    std::string value = std::string(block);
    std::cout << "Should set property " << property << " to value " << value << std::endl;
    
    setProperty(property, value);

    block = white + 1;
    incrementIndent(&block);

    return block;
}

bool Parser::parseNextChunk(char **blockPtr)
{
    char *block = *blockPtr;    

    // just incremented from the first {.
    // we should expect a keyword, or a property next.
    char *space = strchrwhite(block);
    if (space == NULL)
    {
        //std::cout << block << std::endl;
        std::cout << "Nope!" << std::endl;
        return false;
    }

    *space = '\0';
    
    // what are we about to process?

    ParserType type = ParserTypeProperty;
    
    if (strcmp(block, "category") == 0)
    {
        type = ParserTypeObject;
    }
    else if (strcmp(block, "reference") == 0)
    {
        type = ParserTypeReference;
    }
    else if (strcmp(block, "special") == 0)
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

        default:
        break;
    }

    if (block[0] == '}')
    {
        std::cout << "Found a }" << std::endl;
        return false;
    }
   
    return true;
}

void Parser::parse(char *block)
{
    // Get to the beginning of the absolute path...
    incrementIndent(&block);

    char *newline = strchrwhite(block);
    if (newline == NULL) return;

    // We prepare the vessels to accept our data.
    setup();
    
    // Now we can replace the dud path/identifier.
    *newline = '\0';
    std::string path = std::string(block);
    _absolutePath = path;

    std::cout << "Found path " << path << std::endl;
    
    // we also want to isolate the final identifier
    char *slash = strrchr(block, '/');
    slash++;
    std::string identifier = std::string(slash);
    _identifier = identifier;

    std::cout << "Found identifier " << _identifier << std::endl;

    block = newline + 1;
    incrementIndent(&block);
    // we expect the next character to be a {
    if (block[0] != '{')
    {
        std::cout << "Check me: no { for object " << path  << "?" << std::endl;
        return;
    }
    else
    {
        std::cout << "Loading properties for " << _identifier << std::endl;
    }

    std::cout << "Set up ok." << std::endl;
   
    // Get past the {. 
    block++;
    incrementIndent(&block);

    // Loop through the properties.

    bool another = true;

    while (another)
    {    
        another = parseNextChunk(&block);
    }
}

ParserPtr Parser::objectOfType(char *className)
{
    ParserPtr object = ParserPtr();

    if (strcmp(className, "Crystal") == 0) 
    {
        std::cout << "Making object of type " << className << std::endl;
        object = ParserPtr(static_cast<Parser *>(new Crystal()));
    }
    if (strcmp(className, "Polymer") == 0)
    {
        std::cout << "Making object of type " << className << std::endl;
        object = ParserPtr(static_cast<Polymer *>(new Polymer()));
    }
    else
    {
        std::cout << "Do not understand class name " << className << std::endl;
    }

    return object;
}

ParserPtr Parser::processBlock(char *block)
{
    // block should start with class name.
    char *comma = strchr(block, ',');
    *comma = '\0';
    
    ParserPtr object = objectOfType(block);
    if (!object)
    {
        return ParserPtr();
    }

    block = comma + 1;
    object->parse(block);

    return object;
}


