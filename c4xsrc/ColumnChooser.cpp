// cluster4x
// Copyright (C) 2019 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more detailso

#include <iostream>
#include "ColumnChooser.h"
#include "Group.h"
#include "AveVectors.h"
#include <QMessageBox>
#include <algorithm>

ColumnChooser::ColumnChooser()
{
}

void ColumnChooser::addTargetGroup(Group *g)
{
	for (size_t j = 0; j < g->mtzCount(); j++)
	{
		bool found = false;
		for (size_t i = 0; i < _groups.size(); i++)
		{
			for (size_t k = 0; k < _groups[i]->mtzCount(); k++)
			{
				if (g->getMtz(j) == _groups[i]->getMtz(k))
				{
					std::cout << "Warning: Duplicate entry for target!" 
					<< std::endl;
					found = true;
				}
			}
		}
		
		if (!found)
		{
			_mtzs.push_back(g->getMtz(j));
		}
	}

	_groups.push_back(g);
}

double ColumnChooser::evaluate(int enabled, bool work)
{
	double count = 0;
	double sum   = 0;
	
	if (enabled == 1)
	{
		return 100;
	}

	AveVectors *v = Group::topGroup()->getWorkingSet()->vectors;

	for (size_t i = 0; i < _groups.size() - 1; i++)
	{
		for (size_t j = i; j < _groups.size(); j++)
		{
			for (size_t k = 0; k < _groups[i]->mtzCount(); k++)
			{
				MtzFFTPtr fk = _groups[i]->getMtz(k);
				if (_freeMap[fk] == work)
				{
					continue;
				}
				double num = 0;
				double tmp = 0;
				for (size_t l = 0; l < _groups[j]->mtzCount(); l++)
				{
					if (i == j && k == l)
					{
						continue;
					}

					MtzFFTPtr fl = _groups[j]->getMtz(l);
					if (_freeMap[fl] == work)
					{
						continue;
					}

					double add = v->findCorrelation(fk, fl);
					add = std::max(0., add);
					bool same = (i == j);
					tmp += same ? add : 1 - add;
					
					num++;
				}

				tmp /= sqrt(enabled / (enabled - 1));
				
				tmp /= num;
				sum += tmp;
				count++;
			}
		}
	}
	
	return -sum / count;
}

void ColumnChooser::randomFree(double prop)
{
	int expect = _mtzs.size() * prop + 1;
	if (expect < 15)
	{
		expect = 15;
	}
	std::vector<int> shuffled;
	
	for (size_t i = 0; i < _mtzs.size(); i++)
	{
		_freeMap[_mtzs[i]] = false;
		shuffled.push_back(i);
	}

	std::random_shuffle(shuffled.begin(), shuffled.end());

	for (size_t i = 0; i < expect && i < _mtzs.size(); i++)
	{
		_freeMap[_mtzs[shuffled[i]]] = true;
	}
}

int activeColumns()
{
	int enabled = 0;
	for (size_t i = 0; i < AveVectors::titleCount(); i++)
	{
		if (AveVectors::enabled(i))
		{
			enabled++;
		}
	}

	return enabled;
}

void ColumnChooser::enableAllColumns()
{
	for (size_t i = 0; i < AveVectors::titleCount(); i++)
	{
		AveVectors::setEnabled(i, true);
	}
}

bool ColumnChooser::pruneCycle()
{
	AveVectors *v = Group::topGroup()->getWorkingSet()->vectors;
	v->calculate();
	
	int enabled = activeColumns();

	double eval = evaluate(enabled, true);
	
	if (eval > 50)
	{
		return true;
	}

	std::vector<int> remove;
	
	for (size_t i = 0; i < AveVectors::titleCount(); i++)
	{
		if (!AveVectors::enabled(i))
		{
			continue;
		}

		AveVectors::setEnabled(i, false);
		double renewed = evaluate(enabled - 1, true);
		
		if (renewed < eval)
		{
			remove.push_back(i);
		}
		AveVectors::setEnabled(i, true);
	}

	bool changed = false;

	if (enabled - remove.size() > 2)
	{
		for (size_t i = 0; i < remove.size(); i++)
		{
			AveVectors::setEnabled(remove[i], false);
		}
		
		if (remove.size() > 0)
		{
			changed = true;
		}
	}

	return changed;
}

void ColumnChooser::saveColumnSet()
{
	_saved = AveVectors::copyEnabled();
}

void ColumnChooser::loadSavedSet()
{
	AveVectors::setEnabledSet(_saved);
}

void ColumnChooser::prune()
{
	if (_groups.size() == 0)
	{
		std::cout << "No groups marked" << std::endl;
		return;
	}
	
	randomFree(0.5);
	randomFree(0);
	pruneCycle();
	
	double best = 0;
	for (double prop = 0.9; prop > 0.15 && false; prop -= 0.1)
	{
		double swork = 0;
		double sfree = 0;
		double count = 0;
		for (size_t i = 0; i < 12; i++)
		{
			enableAllColumns();
			randomFree(prop);

			bool changed = true;

			do
			{
				changed = pruneCycle();
			}
			while (changed);

			int enabled = activeColumns();
			double work = evaluate(enabled, true);
			double free = evaluate(enabled, false);
			if (work != work || free != free)
			{
				continue;
			}

			swork += work;
			sfree += free;
			count++;
		}
		
		swork /= count;
		sfree /= count;

		std::cout << prop << ", " << swork << ", " << sfree << std::endl;
	}

	QString list;
	list += QString::fromStdString("Most important columns:");
	
	for (size_t i = 0; i < AveVectors::titleCount(); i++)
	{
		if (!AveVectors::enabled(i))
		{
			continue;
		}
		
		list += QString::fromStdString("\n" + AveVectors::title(i));
	}
	
	QMessageBox msg;
	msg.setText(list);
	msg.exec();
}
