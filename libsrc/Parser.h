// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#ifndef __vagabond__Parser__
#define __vagabond__Parser__

#include "shared_ptrs.h"
#include <fstream>

typedef struct
{
    std::string *stringPtr;
    std::string ptrName;
} StringProperty;

class Parser
{
public:
    Parser();

protected:
    virtual std::string getClassName() = 0;
    virtual std::string getIdentifier() = 0;
    virtual void addProperties() = 0;

    void setParent(ParserPtr parent);
    void addStringProperty(std::string className, std::string *ptr);

    void writeToFile(std::ofstream &stream, int indent);
private:
    std::string _className;
    std::string _identifier;
    std::string _absolutePath;    
    ParserPtr _parent;

    std::string getAbsolutePath()
    {
        return _absolutePath;
    }

    std::vector<StringProperty> _stringProperties;

    void makePath();
    void setup();
    bool _setup;
};


#endif
