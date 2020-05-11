//
//  FileReader.h
//  GameDriver
//
//  Created by Helen Ginn on 21/05/2014.
//  Copyright (c) 2014 Helen Ginn. All rights reserved.
//

#ifndef __FileReader__
#define __FileReader__

#include <glob.h> // glob(), globfree()
#include <stdexcept>
#include <cstring>
#include <sstream>

#include <iostream>
#include <string>
#include <vector>
#include <cerrno>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

std::string get_file_contents(std::string filename);

std::vector<std::string> split(const std::string &s, char delim);
bool file_exists(const std::string& name);

std::string getPath(std::string whole);
std::string getFilename(std::string filename);
std::string getBaseFilename(std::string filename);
std::string getBaseFilenameWithPath(std::string filename);
std::string findNextFilename(std::string file);

std::string i_to_str(int val);
std::string f_to_str(double val, int precision);

std::string findNewFolder(std::string prefix = "refine_");

/* Random string things */

void trim(std::string& str);
void to_lower(std::string &str);
void to_upper(std::string &str);

void print_cc_diff(double diff, int limit);

inline std::vector<std::string> glob(const std::string& pattern) {
	using namespace std;

	// glob struct resides on the stack
	glob_t glob_result;
	memset(&glob_result, 0, sizeof(glob_result));

	// do the glob operation
	int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
	if (return_value != 0) 
	{
		globfree(&glob_result);

		if (return_value == GLOB_NOMATCH)
		{
			return std::vector<std::string>();
		}

		stringstream ss;
		ss << "glob() failed with return_value " << return_value << endl;
		throw std::runtime_error(ss.str());
	}

	// collect all the filenames into a std::list<std::string>
	vector<string> filenames;

	for (size_t i = 0; i < glob_result.gl_pathc; i++)
	{
		filenames.push_back(string(glob_result.gl_pathv[i]));
	}

	// cleanup
	globfree(&glob_result);

	// done
	return filenames;
}

class FileReader
{

public:
	static void makeDirectoryIfNeeded(std::string _dir)
	{
		DIR *dir = opendir(_dir.c_str());

		if (dir)
		{
			closedir(dir);
		}
		else if (ENOENT == errno)
		{
			mkdir(_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		}
	}

	static void setOutputDirectory(std::string _dir)
	{
		outputDir = _dir;

		makeDirectoryIfNeeded(outputDir);
	}

	static std::string addOutputDirectory(std::string filename,
	                                      std::string subdir = "")
	{
		if (!outputDir.length())
		{
			return filename;
		}

		if (outputDir[0] == '/')
		{
			return outputDir + "/" + filename;
		}
		
		if (subdir.length())
		{
			makeDirectoryIfNeeded("./" + outputDir + "/" + subdir);
			subdir += "/";
		}

		std::string fullPath = "./" + outputDir + "/" + subdir + filename;
		return fullPath;
	}

private:
	static std::string outputDir;

};

#endif /* defined(__GameDriver__FileReader__) */
