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

std::vector<std::string> &FileReader::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> FileReader::split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

bool FileReader::exists(const std::string& name)
{
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}

std::string FileReader::get_file_contents(const char *filename)
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

void read_file(const char *name, char **buffer, int *length)
{
    FILE *file;
    unsigned long file_length;
    
    file = fopen(name, "rb");
    if (!file)
    {
        fprintf(stderr, "Unable to open file %s", name);
        return;
    }
    
    /* Get file length, return to beginning */
    fseek(file, 0, SEEK_END);
    file_length = ftell(file);
    *length = (int)file_length;
    fseek(file, 0, SEEK_SET);
    
    /* Allocate memory */
    (*buffer) = (char *) malloc(file_length + 1);
    
    if (!(*buffer))
    {
        fprintf(stderr, "Memory error!");
        fclose(file);
        return;
    }
    
    /* Read file contents into buffer */
    fread(*buffer, file_length, 1, file);
    fclose(file);
}
