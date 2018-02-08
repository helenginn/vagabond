// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#include "Parser.h"
#include <sstream>
#include <iostream>
#include <iomanip>

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
    child->setParent(this);
    _parserList[category].push_back(child);
    child->setup();
}

std::string indent(int num)
{
    std::ostringstream stream;

    for (int i = 0; i < num; i++)
    {
        stream << "  ";
    }
    return stream.str();
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
        stream << indent(in) << name << " = \"" << *ptr << "\"" << std::endl;
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
        if (!ptr || delegate) continue;

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
        int *ptr = _boolProperties[i].boolPtr;
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

void Parser::writeToFile(std::ofstream &stream, int in)
{
    setup();

    stream << "vagabond data structure v0.0" << std::endl;

    outputContents(stream, in);
}


void Parser::addReference(std::string category, ParserPtr cousin)
{
    _referenceList[category].push_back(cousin);
    
}

