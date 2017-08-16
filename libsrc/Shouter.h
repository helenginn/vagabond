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

void shout_at_user(std::string fix_me_message);
void shout_at_helen(std::string fix_me_message);
void warn_user(std::string cautionary_tale);
void shout_timer(time_t wall_start, std::string job);

#endif /* defined(__vagabond__Shouter__) */
