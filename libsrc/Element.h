//
//  Element.h
//  vagabond
//
//  Created by Helen Ginn on 16/07/2017.
//  Copyright (c) 2017 Strubi. All rights reserved.
//

#ifndef __vagabond__Element__
#define __vagabond__Element__

#include <stdio.h>
#include <string>
#include <vector>
#include "shared_ptrs.h"
#include "Distributor.h"

class Element : public Distributor
{
public:
    static void setupElements();

    Element(std::string symbol, std::string name, double electrons,
            const float *scatter);
    static ElementPtr getElement(std::string symbol);
    virtual FFTPtr getDistribution(bool quick = false);

    std::string *getSymbolPtr()
    {
        return &_symbol;
    }

    std::string getSymbol()
    {
        return _symbol;
    }

    std::string getName()
    {
        return _name;
    }

    double electronCount()
    {
        return _electrons;
    }
    
protected:
    static double getVoxelValue(void *obj, double x, double y, double z);

private:
    std::string _symbol;
    std::string _name;
    double _electrons;
    
    float _scattering[62];
    FFTPtr _fft;

    static std::vector<ElementPtr> elements;

    void setSymbol(std::string symbol)
    {
        _symbol = symbol;
    }

    void setName(std::string name)
    {
        _name = name;
    }

};

#endif /* defined(__vagabond__Element__) */
