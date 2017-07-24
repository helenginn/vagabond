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

std::string get_file_contents(std::string filename);

std::vector<std::string> split(const std::string &s, char delim);
bool file_exists(const std::string& name);

std::string getFilename(std::string filename);
std::string getBaseFilename(std::string filename);

/* Random string things */

void trim(std::string& str);
void to_lower(std::string &str);

#endif /* defined(__GameDriver__FileReader__) */
