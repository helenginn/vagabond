//
//  Shouter.cpp
//  vagabond
//
//  Created by Helen Ginn on 13/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Shouter.h"
#include <iostream>
#include <time.h>

void shout_at_user(std::string fix_me_message)
{
	std::cout << "****************************************" << std::endl;
	std::cout << "**               Error                **" << std::endl;
	std::cout << "****************************************" << std::endl;
	std::cout << std::endl;
	std::cout << fix_me_message << std::endl;
	std::cout << std::endl;
	std::cout << "****************************************" << std::endl;
	std::cout << "**      Please fix me and re-run      **" << std::endl;
	std::cout << "****************************************" << std::endl;

	exit(1);
}

Shouter::Shouter(std::string message)
{
	_message = message;
}

void Shouter::shoutToStdOut()
{
	shout_at_user(_message);
}

void shout_at_helen(std::string fix_me_message)
{
	std::cout << "****************************************" << std::endl;
	std::cout << "**           Bug detected!            **" << std::endl;
	std::cout << "****************************************" << std::endl;
	std::cout << std::endl;
	std::cout << fix_me_message << std::endl;
	std::cout << std::endl;
	std::cout << "****************************************" << std::endl;
	std::cout << "**     Please contact the author      **" << std::endl;
	std::cout << "**     and provide this log file      **" << std::endl;
	std::cout << "****************************************" << std::endl;

	exit(1);
}

void warn_user(std::string cautionary_tale)
{
	std::cout << "*** Warning: " << cautionary_tale << " ***\n" << std::endl;
}

void shout_timer(time_t wall_start, std::string job)
{
	time_t wall_end;
	time(&wall_end);

	time_t diff = wall_end - wall_start;
	time_t seconds = (diff % 60);
	time_t minutes = (diff - seconds) / 60;

	std::cout << "~ Clock time for " << job << ": " << minutes 
	<< "m" << seconds << "s. ~" << std::endl;
	std::cout << std::endl;
}
