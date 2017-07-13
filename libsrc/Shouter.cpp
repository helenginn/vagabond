//
//  Shouter.cpp
//  vagabond
//
//  Created by Helen Ginn on 13/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#include "Shouter.h"
#include <iostream>

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