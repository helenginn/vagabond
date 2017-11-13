//
//  FileReader.cpp
//  GameDriver
//
//  Created by Helen Ginn on 21/05/2014.
//  Copyright (c) 2014 Helen Ginn. All rights reserved.
//

#include "FileReader.h"
#include <fstream>
#include <sstream>
#include <cerrno>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <iomanip>
#include <algorithm>

std::string FileReader::outputDir;

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

bool file_exists(const std::string& name)
{
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}

std::string get_file_contents(std::string filename)
{
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    
    if (in)
    {
        std::string contents;
        in.seekg(0, std::ios::end);
        contents.resize((unsigned long)in.tellg());
        in.seekg(0, std::ios::beg);
        in.read(&contents[0], contents.size());
        in.close();
        return(contents);
    }
	

    std::string errString = "Could not get file contents for file " + std::string(filename);
	std::cout << errString << std::endl;

    throw(errno);
}


std::string getFilename(std::string filename)
{
	size_t pos = filename.rfind("/");
	if(pos == std::string::npos)  //No path.
		return filename;

	return filename.substr(pos + 1, filename.length());
}

std::string getBaseFilename(std::string filename)
{
	std::string fName = getFilename(filename);
	size_t pos = fName.rfind(".");
	if(pos == std::string::npos)  //No extension.
		return fName;

	if(pos == 0)    //. is at the front. Not an extension.
		return fName;

	return fName.substr(0, pos);
}

void trim(std::string &str)
{
	std::string::size_type pos = str.find_last_not_of(' ');
	if(pos != std::string::npos)
	{
		str.erase(pos + 1);
		pos = str.find_first_not_of(' ');
		if(pos != std::string::npos) str.erase(0, pos);
	}
	else str.erase(str.begin(), str.end());
}

void to_lower(std::string &str)
{
	std::transform(str.begin(), str.end(), str.begin(), ::tolower);
}

void to_upper(std::string &str)
{
	std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

std::string i_to_str(int val)
{
	std::ostringstream ss;
	ss << val;
	std::string temp = ss.str();

	return temp;
}

std::string f_to_str(double val, int precision)
{
	std::ostringstream ss;
	if (precision > 0)
	{
		ss << std::fixed << std::setprecision(precision);
	}
	else if (precision < 0)
	{
		ss << std::fixed;
	}

	ss << val;
	std::string temp = ss.str();

	return temp;
}
