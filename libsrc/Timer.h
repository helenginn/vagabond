//
//  Timer.h
//  vagabond
//
//  Created by Helen Ginn on 13/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Timer__
#define __vagabond__Timer__

#include <time.h>
#include <string>

/**
 * \class Timer
 * \brief Start, stop and report timing in minutes/seconds within other
 * areas of Vagabond.
 */

class Timer
{
public:
	Timer();	
	Timer(std::string name, bool start = false);	
	
	void start();
	void stop(bool quick = false);

	
	void quickReport();
	void report();
private:
	time_t wall;
	time_t accumulative;
	std::string _name;
};


#endif

