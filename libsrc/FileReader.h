//
//  FileReader.h
//  GameDriver
//
//  Created by Helen Ginn on 21/05/2014.
//  Copyright (c) 2014 Helen Ginn. All rights reserved.
//

#ifndef __FileReader__
#define __FileReader__

#include <iostream>
#include <string>
#include <vector>

namespace FileReader
{
    std::string get_file_contents(const char *filename);
    
    std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
    std::vector<std::string> split(const std::string &s, char delim);
    bool exists(const std::string& name);
};

void read_file(const char *name, char **buffer, int *length);
std::string getFilename(std::string filename);
std::string getBaseFilename(std::string filename);

#endif /* defined(__GameDriver__FileReader__) */
