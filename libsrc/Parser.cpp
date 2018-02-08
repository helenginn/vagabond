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
}

void Parser::setup()
{
    if (_setup) return;

    _identifier = getIdentifier(); 
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

void Parser::setParent(ParserPtr parent)
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

std::string indent(int num)
{
    std::ostringstream stream;

    for (int i = 0; i < num; i++)
    {
        stream << "  ";
    }
    return stream.str();
}

void Parser::writeToFile(std::ofstream &stream, int i)
{
    setup();
    stream << std::setprecision(5);
    stream << indent(i) << "object " << _className
           << ", " << _absolutePath << std::endl;
    stream << indent(i) << "{" << std::endl;
    i++;
    
    for (int i = 0; i < _stringProperties.size(); i++)
    {
        std::string name = _stringProperties[i].ptrName;
        std::string *ptr = _stringProperties[i].stringPtr;
        if (!ptr) continue;
        stream << indent(i) << name << " = \"" << *ptr << "\"" << std::endl;
    }
 
    for (int i = 0; i < _doubleProperties.size(); i++)
    {
        std::string name = _doubleProperties[i].ptrName;
        double *ptr = _doubleProperties[i].doublePtr;
        if (!ptr) continue;
        stream << indent(i) << name << " = " << *ptr << std::endl;
    }

    stream << indent(i) << "stuff..." << std::endl;
    i--;
    stream << indent(i) << "}" << std::endl;
}
