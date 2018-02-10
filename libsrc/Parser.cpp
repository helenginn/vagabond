// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#include "Parser.h"
#include <iostream>
#include <iomanip>
#include "Crystal.h"
#include "Polymer.h"
#include "Atom.h"
#include "Bond.h"
#include "Monomer.h"
#include "Absolute.h"
#include "Sidechain.h"
#include "Backbone.h"

ParserMap Parser::_allParsers;

Parser::Parser()
{
    _setup = false;
    _parent = NULL;
}

void Parser::setup(bool isNew)
{
    if (_setup) return;
    
    if (!isNew)
    {
        _identifier = getParserIdentifier(); 
        _className = getClassName();

        makePath();
    }
    
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
                               void *delegate, Encoder encoder,
                               Decoder decoder) 
{
    CustomProperty property;
    property.ptrName = className;
    property.objPtr = ptr;
    property.delegate = delegate;
    property.encoder = encoder;
    property.decoder = decoder;
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
        stream << indent(in) << name << " = " << *ptr << std::endl;
    }

    for (int i = 0; i < _vec3Properties.size(); i++)
    {
        std::string name = _vec3Properties[i].ptrName;
        vec3 *ptr = _vec3Properties[i].vec3Ptr;
        if (!ptr) continue;
        stream << indent(in) << name << " = " << ptr->x
               << "," << ptr->y << "," << ptr->z << std::endl;
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
    _allParsers.clear();
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

char *Parser::parseNextSpecial(char *block)
{
    // posied at the word just after 'special'.

    char *white = strchrwhite(block);
    *white = 0;

    std::string specialName = std::string(block);
//    std::cout << "Special name is " << specialName << std::endl;
    block = white + 1;
    incrementIndent(&block);

    if (block[0] != '{')
    {
        std::cout << "Something's wrong - was expecting a {!" << std::endl;
        return white;
    }

    block++;
    incrementIndent(&block);
    
    for (int i = 0; i < _customProperties.size(); i++)
    {
        CustomProperty property = _customProperties[i];
        if (property.ptrName == specialName)
        {
            Decoder decoder = property.decoder;
            void *delegate = property.delegate;
            void *ptr = property.objPtr;
            block = (*decoder)(delegate, ptr, block);
        }
    }
 
    if (block[0] != '}')
    {
        std::cout << "Why has this special thing not ended?" << std::endl;
        std::cout << "Actual char: " << block[0] << std::endl;
    }

    block++;
    incrementIndent(&block);

    return block;
}

char *Parser::parseNextReference(char *block)
{
    // poised at the word just after 'references'.

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

    // expecting a fuckton of references now
    // will go through later and fill them in

    while (true)
    {
        white = strchrwhite(block);
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
    
    if (strcmp(block, "object") != 0)
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
        _allParsers[path] = object;

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
        
        if (strcmp(block, "object") != 0)
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
//        std::cout << "Found a } after end of category" << std::endl;
        block++;
        incrementIndent(&block);
    }

    return block;
}

void Parser::setProperty(std::string property, std::string value)
{
//    std::cout << "Setting property " << property << " to " << value << std::endl;

    for (int i = 0; i < _stringProperties.size(); i++)
    {
        if (_stringProperties[i].ptrName == property)
        {
            if (value == "__NULL__")
            {
                *_stringProperties[i].stringPtr = "";
                return; 
            }

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
        std::cout << "Nope!" << std::endl;
        std::cout << block << std::endl;
        return false;
    }

    *space = '\0';
    
    // what are we about to process?

    ParserType type = ParserTypeProperty;
    
    if (strcmp(block, "category") == 0)
    {
        type = ParserTypeObject;
    }
    else if (strcmp(block, "references") == 0)
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

        case ParserTypeObject:
        *blockPtr = parseNextObject(block);
        break;

        case ParserTypeSpecial:
        *blockPtr = parseNextSpecial(block);
        break;

        case ParserTypeReference:
        *blockPtr = parseNextReference(block);
        break;

        default:
        break;
    }

    block = *blockPtr;

//    std::cout << "Completed a 'next' thing, char '" << block[0] << "'" << std::endl;

    if (*blockPtr == NULL)
    {
        std::cout << "Parsing error occurred." << std::endl;
        return false;
    }

    if (block[0] == '}')
    {
//        std::cout << "Found a }" << std::endl;
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

//    std::cout << "Found path " << path << std::endl;
    
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
    else
    {
//        std::cout << "Loading properties for " << _identifier << std::endl;
    }

//    std::cout << "Set up ok." << std::endl;
   
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
        object = ParserPtr(static_cast<Parser *>(new Crystal()));
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
    else if (strcmp(className, "Bond") == 0)
    {
        object = ParserPtr(static_cast<Bond *>(new Bond()));        
    }
    else if (strcmp(className, "Molecule") == 0)
    {
        object = ParserPtr(static_cast<Molecule *>(new Molecule()));        
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

ParserPtr Parser::processBlock(char *block)
{
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

    // Resolve dangling references.
    object->resolveReferences();

    // Loop through all objects to allow them to finish up.
    for (ParserMap::iterator it = _allParsers.begin();
         it != _allParsers.end(); it++)
    {
        ParserPtr aParser = it->second;
        aParser->postParseTidy();
    }

    if (success != NULL)
    {
        return object;
    }
    
    return ParserPtr();
}

ParserPtr Parser::resolveReference(std::string reference)
{
    return _allParsers[reference];
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
            it->second[j]->resolveReferences();
        }
    }
}

