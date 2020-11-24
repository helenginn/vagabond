//
//  Timer.cpp
//  vagabond
//
//  Created by Helen Ginn on 13/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Timer.h"
#include <iostream>
#include <iomanip>

Timer::Timer()
{
	wall = 0;
	accumulative = 0;
	_milliseconds = 0;
	_fine = false;
	start();	
}

Timer::Timer(std::string name, bool _start)
{
	_name = name;
	wall = 0;
	accumulative = 0;
	_fine = false;
	_milliseconds = 0;
	
	if (_start)
	{
		start();	
	}
}


void Timer::start()
{
	time(&wall);
	_tStart = std::chrono::high_resolution_clock::now();
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

	_tEnd = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> time_span = _tEnd - _tStart;
	_milliseconds += time_span.count();
}

void Timer::report()
{
	if (wall > 0)
	{
		stop();
	}

	if (_fine)
	{
		std::cout << "~ Finely measured time for " << _name << ": ";
		std::cout <<  std::setprecision(3) << _milliseconds / 1000;
		std::cout << "s. ~" << std::endl;
		std::cout << std::endl;
	}
	else
	{
		time_t seconds = (accumulative % 60);
		time_t minutes = (accumulative - seconds) / 60;

		std::cout << "~ Clock time for " << _name << ": " << minutes << "m" << seconds << "s. ~" << std::endl;
		std::cout << std::endl;
	}
}
