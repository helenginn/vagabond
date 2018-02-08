// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#include "Parser.h"
#include <sstream>
#include <iostream>

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
    stream << indent(i) << "object " << _className
           << ", " << _absolutePath << std::endl;
    stream << indent(i) << "{" << std::endl;
    i++;
    stream << indent(i) << "stuff..." << std::endl;
    i--;
    stream << indent(i) << "}" << std::endl;
}
