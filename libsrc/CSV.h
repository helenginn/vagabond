//
//  CSV.h
//   cppxfel - a collection of processing algorithms for XFEL diffraction data.

//    Copyright (C) 2017  Helen Ginn
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.


#ifndef __cppxfel__CSV__
#define __cppxfel__CSV__

#include <stdio.h>
#include <vector>
#include <string>
#include "shared_ptrs.h"
#include <cstdarg>
#include <map>

typedef enum
{
    PlotVerticalLine = '|',
    PlotHorizontalLine = '_',
    PlotHorizontalTickMark = '-',
    PlotVerticalTickMark = '\'',
    PlotBlank = ' ',
    
} PlotChar;

typedef enum
{
    ConvolutionTypeSuperGaussian,
    ConvolutionTypeUniform,
} ConvolutionType;

typedef std::vector<double> Entry;
typedef std::map<int, char> Row;
typedef std::map<int, Row > Plot;

class CSV
{
private:
    std::vector<std::string> headers;
    std::vector<Entry> entries;
    std::vector<double> histCategories;
    double minX, minY, maxX, maxY;
    bool didSetMinMaxXY;
    
    std::string mapToAscii(Plot plot);
    void writeStringToPlot(std::string text, Plot *plot, int x, int y);
public:
    CSV()
    {
        
    }

    CSV(int count, ...)
    {
        didSetMinMaxXY = false;
        
        va_list arguments;
        va_start(arguments, count);
        
        for (int i = 0; i < count; i++)
        {
            std::string header = std::string(va_arg(arguments, char *));
            addHeader(header);
        }
        
        va_end(arguments);
    }


	void setupHistogram(double start, double end, double interval, std::string catHeader, int count, ...)
	{
		va_list arguments;
		va_start(arguments, count);

		headers.push_back(catHeader);

		for (int i = 0; i < count; i++)
		{
			std::string header = std::string(va_arg(arguments, char *));
			headers.push_back(header);
		}

		va_end(arguments);

		double value = start;
		while (value <= end)
		{
			addPartialEntry(1, value);
			histCategories.push_back(value);
			value += interval;
		}
	}

	void minMaxCol(int col, double *min, double *max, bool round = false);
    void addOneToFrequency(double category, std::string whichHeader, double weight = 1, std::string categoryHeader = "");
    void addOneToFrequency(double category, int column, double weight = 1, int categoryNum = 0);
    int findHeader(std::string whichHeader);

    ~CSV();
    
    void plotPNG(std::map<std::string, std::string> properties);
	
    void addPartialEntry(int dummy, ...);
    void addEntry(int dummy, ...);
    void writeToFile(std::string filename);
    double valueForEntry(std::string header, int entry);
    double valueForHistogramEntry(std::string whichHeader, double value, std::string categoryHeader = "");
    double valueForHistogramEntry(int whichHeader, double value, int categoryHeader = 0);
    void histogram(std::map<double, int> histogram);
    std::string plotColumns(int col1, int col2);
    void resetColumn(std::string header, double value = 0);
    
    void setValueForEntry(int entry, std::string header, double value);

	void reserveEntries(unsigned long num)
	{
		entries.reserve(num);
	}

	void addEntry(std::vector<double> entry)
    {
        entries.push_back(entry);
    }
    
    int entryCount()
    {
        return (int)entries.size();
    }

	std::vector<double> entry(int i)
	{
		return entries[i];
	}

	double valueForEntry(int headerNum, int entry)
	{
		return entries[entry][headerNum];
	}

    int headerCount()
    {
        return (int)headers.size();
    }
    
    void addHeader(std::string header)
    {
        headers.push_back(header);
    }
    
    void setMinMaxXY(double _minX, double _minY, double _maxX, double _maxY)
    {
        minX = _minX;
        minY = _minY;
        maxX = _maxX;
        maxY = _maxY;
        didSetMinMaxXY = true;
    }
};

#endif /* defined(__cppxfel__CSV__) */
