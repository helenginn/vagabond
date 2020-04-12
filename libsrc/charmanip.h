// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#ifndef __vagabond__charmanip__
#define __vagabond__charmanip__

#include <sstream>
#include <cstring>
#include <iostream>

inline std::string indent(int num)
{
	std::ostringstream stream;

	for (int i = 0; i < num; i++)
	{
		stream << "  ";
	}
	return stream.str();
}

inline char *strchrwhite(char *block)
{
	while (true)
	{
		if (*block == ' ' || *block == '\n' 
		    || *block == '\t' || *block == '\0')
		{
			return block;
		}
		
		block++;
	}
	
	return NULL;
}

inline void incrementIndent(char **block)
{
	while ((*block)[0] == ' ' || (*block)[0] == '\t' || (*block)[0] == '\n' || (*block)[0] == '\0')
	{
		(*block)++;
	}
}

inline char *keywordValue(char *block, char **keyword, char **value) 
{
	char *space = strchrwhite(block);

	if (space == NULL)
	{
		std::cout << "Space is just null" << std::endl;
		return NULL;
	}

	*space = '\0';
	*keyword = block;
	block = space + 1;
	incrementIndent(&block);

	// Don't panic, we probably just have an 'object'.
	if (block[0] != '=')
	{
		//        std::cout << "keyword: " << *keyword << " - block char " << *block << std::endl;
		return block;
	}

	block++;
	incrementIndent(&block);

	space = strchrwhite(block);
	*space = '\0';

	*value = block;
	block = space + 1;
	incrementIndent(&block);

	return block;
}

#endif
