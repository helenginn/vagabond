// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#include "VBondReader.h"
#include "FileReader.h"
#include "Crystal.h"

CrystalPtr VBondReader::getCrystal()
{
    std::string vbondStr = get_file_contents(_filename);
  
    char *vbond = &vbondStr[0];
    char header[] = "vagabond data structure";
    int headLength = strlen(header);
    vbond[headLength] = '\0';

    if (strcmp(vbond, header) != 0)
    {
        std::cout << "Aborting: " << _filename << " does not look like a vagabond data structure." << std::endl;
        return CrystalPtr();
    }

    vbond = &vbond[headLength + 1];
    char *endline = strchr(vbond, '\n');
    *endline = '\0'; endline++;
    std::cout << "Vagabond data structure file is version " << vbond << std::endl;
    vbond = endline;
    char *space = strchr(vbond, ' ');
    *space = '\0';   
 
    if (strcmp(vbond, "object") == 0)
    {
        vbond = space + 1;
        ParserPtr parser = Parser::processBlock(vbond);

        if (!parser)
        {
            std::cout << "Parser could not make top level object." << std::endl;
            return CrystalPtr();
        }

        std::string className = parser->getClassName();
        
        if (className == "Crystal")
        {
            CrystalPtr crystal = ToCrystalPtr(parser);
            return crystal;
        }
        else
        {
            std::cout << "Was expecting a Crystal object." << std::endl;
            return CrystalPtr();
        }
    }
    else
    {
        std::cout << "Was expecting object in " << _filename << std::endl;
        std::cout << "Instead, got " << vbond << std::endl;
    }
 
    return CrystalPtr();
}

VBondReader::VBondReader()
{

}
