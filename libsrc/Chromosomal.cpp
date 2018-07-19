//
//  Chromosomal.cpp
//  vagabond
//
//  Created by Helen Ginn on 01/07/2018.
//  Copyright (c) 2018 Helen Ginn. All rights reserved.
//

#include <iostream>
#include "Chromosomal.h"
#include <map>
#include <algorithm>

typedef struct
{
	ChromosomalPtr specimen;
	double score;
} ScorePair;

bool greater_score(ScorePair &a, ScorePair &b)
{
	return a.score > b.score;
}

Chromosomal::Chromosomal()
{
	_sampled = 50;
}

void Chromosomal::testPopulation()
{
	std::vector<ScorePair> scores;	

	for (int i = 0; i < _copies.size(); i++)	
	{
		double result = _copies[i]->evaluate();
		ScorePair score;
		score.specimen = _copies[i];
		score.score = -result;
		scores.push_back(score);
		std::cout << "." << std::flush;
	}

	std::cout << std::endl;
	
	std::sort(scores.begin(), scores.end(), greater_score);
	_survivors.clear();

	for (int i = 0; i < scores.size() && i < _sampled / 10; i++)	
	{
		double result = scores[i].score;
		ChromosomalPtr specimen = scores[i].specimen;
		_survivors.push_back(specimen);
		std::cout << "Top: " << result << " ";
		specimen->geneticCode();
	}

	_survivors[0]->evaluate();
	_copies.clear();
}

void Chromosomal::breed()
{
	while (_copies.size() < _sampled)
	{
		for (int i = 0; i < _survivors.size() - 1; i++)
		{
			ChromosomalPtr copy = _survivors[i]->makeCopy();
			_copies.push_back(copy);

			for (int j = i; j < _survivors.size(); j++)
			{
				ChromosomalPtr copy = _survivors[i]->makeCopy();
				copy->mutate();
				_copies.push_back(copy);
			}
		}
	}
	
	for (int i = 0; i < _sampled / 4; i++)
	{
		int r1 = rand() % _copies.size();
		int r2 = rand() % _copies.size();
		
		if (r1 == r2)
		{
			i--;
			continue;
		}
		
		_copies[r1]->haveSexWith(&*_copies[r2]);
	}
	
	std::cout << "Next lot: " << _copies.size() << std::endl;
}

void Chromosomal::evolve()
{
	geneticCode(); 
	for (int i = 0; i < _sampled; i++)
	{
		ChromosomalPtr copy = makeCopy();
		copy->randomise();
		_copies.push_back(copy);
	}

	for (int i = 0; i < 1; i++)	
	{
		testPopulation();
		breed();
	}

	testPopulation();
	
	copyOver(&*_survivors[0]);
	geneticCode();
}
