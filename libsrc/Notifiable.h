//
//  Notifiable.h
//  Notifiable
//
//  Created by Helen Ginn on 21/01/2018
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#ifndef __Vagabond__Notifiable__
#define __Vagabond__Notifiable__

#include <string>

class Notifiable
{
public:
    Notifiable()
    {
        _enabled = false;
    }
    
    virtual void enable()
    {
        _enabled = true;
    }
    
    virtual void disable()
    {
        _enabled = false;
    }

    virtual void setMessage(std::string message)
    {
        _message = message;
    }

private:
    std::string _message;
    bool _enabled;
};

#endif
