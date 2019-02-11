//
//  Timer.cpp
//  vagabond
//
//  Created by Helen Ginn on 13/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Timer.h"
#include <iostream>

Timer::Timer()
{
	wall = 0;
	accumulative = 0;
	start();	
}

Timer::Timer(std::string name, bool _start)
{
	_name = name;
	wall = 0;
	accumulative = 0;
	
	if (_start)
	{
		start();	
	}
}


void Timer::start()
{
	time(&wall);
}

void Timer::quickReport()
{
	if (wall > 0)
	{
		stop();
	}
	
	time_t seconds = (accumulative % 60);
	time_t minutes = (accumulative - seconds) / 60;

	std::cout << "(";
	
	if (minutes > 0)
	{
		std::cout << minutes << "m ";
	}
	
	std::cout << seconds << "s)";
}

void Timer::stop(bool quick)
{
	time_t wall_end;
	time(&wall_end);

	time_t diff = wall_end - wall;

	wall = 0;
	
	accumulative += diff;
}

void Timer::report()
{
	if (wall > 0)
	{
		stop();
	}
	
	time_t seconds = (accumulative % 60);
	time_t minutes = (accumulative - seconds) / 60;

	std::cout << "~ Clock time for " << _name << ": " << minutes << "m" << seconds << "s. ~" << std::endl;
	std::cout << std::endl;
}
