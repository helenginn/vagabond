//
//  Shouter.h
//  vagabond
//
//  Created by Helen Ginn on 13/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Shouter__
#define __vagabond__Shouter__

#include <stdio.h>
#include <string>
#include <vector>

void shout_at_user(std::string fix_me_message);
void shout_at_helen(std::string fix_me_message);
void warn_user(std::string cautionary_tale);
void shout_timer(time_t wall_start, std::string job);

typedef struct
{
	std::vector<std::string> availabels;
	std::vector<std::string> types;
	std::string wanted;
	std::string original;
} LabelChoice;

class Shouter
{
public:
	Shouter(std::string message);

	std::string getMessage()
	{
		return _message;
	}
	
	void setChoice(LabelChoice choice)
	{
		_choice = choice;
		_fixable = true;
	}
	
	bool isFixable()
	{
		return _fixable;
	}
	
	LabelChoice &getChoice()
	{
		return _choice;
	}
	
	static void throw_shout(std::string message);
	void shoutToStdOut();

private:
	LabelChoice _choice;
	
	bool _fixable;
	std::string _message;
};

#endif /* defined(__vagabond__Shouter__) */
