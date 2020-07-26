// Vagabond
// Copyright (C) 2017-2018 Helen Ginn
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
// Please email: vagabond @ hginn.co.uk for more details.


#include "CSV.h"
#include <fstream>
#include <float.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <algorithm>
#include "FileReader.h"
#include "PNGFile.h"
#include <cmath>
#include "vec3.h"

int CSV::findHeader(std::string whichHeader)
{
	int chosenHeader = -1;

	for (int i = 0; i < headers.size(); i++)
	{
		if (headers[i] == whichHeader)
		{
			chosenHeader = i;
			break;
		}
	}


	if (chosenHeader == -1)
	{
		std::cout << "Error: " << whichHeader << " does not exist in table." << std::endl;
	}

	return chosenHeader;
}

void CSV::addPartialEntry(int dummy, ...)
{
	va_list arguments;
	va_start(arguments, dummy);
	Entry newEntry;

	for (int i = 0; i < dummy; i++)
	{
		double value = va_arg(arguments, double);
		newEntry.push_back(value);
	}

	int difference = (int)(headers.size() - newEntry.size());

	while (difference > 0)
	{
		newEntry.push_back(0);
		difference--;
	}

	entries.push_back(newEntry);
}

void CSV::addEntry(int dummy, ...)
{
	va_list arguments;
	va_start(arguments, dummy);
	Entry newEntry;

	for (int i = 0; i < headers.size(); i++)
	{
		double value = va_arg(arguments, double);
		newEntry.push_back(value);
	}

	entries.push_back(newEntry);
}

void CSV::writeToFile(std::string filename)
{
	std::ofstream csv;
	std::string outputFile = FileReader::addOutputDirectory(filename, _subdir);

	csv.open(outputFile.c_str());

	for (int i = 0; i < headers.size(); i++)
	{
		csv << headers[i] << ", ";
	}

	csv << std::endl;

	for (int i = 0; i < entries.size(); i++)
	{
		Entry anEntry = entries[i];

		for (int j = 0; j < anEntry.size(); j++)
		{
			csv << std::setprecision(14) << anEntry[j] << ", ";
		}

		csv << std::endl;
	}

	csv.close();
}

double CSV::valueForEntry(std::string header, int entry)
{
	for (int i = 0; i < headerCount(); i++)
	{
		if (headers[i] == header)
		{
			return entries[entry][i];
		}
	}

	return 0;
}

std::string CSV::mapToAscii(Plot plot)
{
	std::ostringstream stream;

	for (Plot::iterator it = plot.begin(); it != plot.end(); it++)
	{
		Row row = it->second;

		stream << "N: ";

		for (Row::iterator it2 = row.begin(); it2 != row.end(); it2++)
		{
			char charray[2];
			charray[0] = it2->second;
			charray[1] = '\0';

			stream << charray;
		}

		stream << std::endl;
	}

	return stream.str();
}

void CSV::minMaxCol(int col, double *min, double *max, bool round)
{
	*min = FLT_MAX;
	*max = -FLT_MAX;

	for (int i = 0; i < entries.size(); i++)
	{
		if (entries[i][col] > *max)
		*max = entries[i][col];

		if (entries[i][col] < *min)
		*min = entries[i][col];
	}

	if (round)
	{
		double logMax = log10(*max);

		double whole = floor(logMax);
		double remainder = fmod(logMax, 1);

		for (int i = 1; i <= 10; i++)
		{
			const double logi = log10(i);
			if (remainder <= logi)
			{
				whole += logi;
				break;
			}
		}

		*max = pow(10, whole);
	}
}

void CSV::writeStringToPlot(std::string text, Plot *plot, int y, int x)
{
	if (!text.length())
	{
		return;
	}

	for (int i = 0; i < text.length(); i++)
	{
		(*plot)[y][x + i] = text[i];
	}
}

typedef enum
{
	GraphStyleScatter,
	GraphStyleLine,
	GraphStyleHeatMap,
} GraphStyle;

void CSV::plotPNG(std::map<std::string, std::string> properties)
{
	int count = 0;

	double height = 900.;

	if (properties.count("height"))
	{
		height = atof(properties["height"].c_str());
	}

	double width = 1200.;

	if (properties.count("width"))
	{
		width = atof(properties["width"].c_str());
	}

	const int xIntervals = 5;
	const int yIntervals = 5;
	const double xAxis = 0.8;
	const double yAxis = 0.7;
	double fastStride = 0;

	if (properties.count("stride"))
	{
		fastStride = atoi(properties["stride"].c_str());
	}

	std::string filename = "graph_test.png";

	if (properties.count("filename"))
	{
		filename = properties["filename"] + ".png";
	}

	PNGFilePtr png = PNGFilePtr(new PNGFile(filename, width, height));
	png->setSubDirectory(_subdir);
	png->setPlain();
	png->setCentre(width * (1 - xAxis) / 2, height * (yAxis + 1) / 2);
	png->drawLine(0, 0, width * xAxis, 0, 0, 0, 0, 0);
	png->drawLine(0, 0, 0, -height * yAxis, 0, 0, 0, 0);

	while (true)
	{
		std::string headerKey = "Header" + i_to_str(count);
		if (!properties.count("x" + headerKey) || !properties.count("y" + headerKey))
		{
			break;
		}

		std::string xHeader = properties["x" + headerKey];
		std::string yHeader = properties["y" + headerKey];

		int xCol = findHeader(xHeader);
		int yCol = findHeader(yHeader);

		int zCol = 0;
		bool hasZ = properties.count("z" + headerKey);

		if (hasZ)
		{
			std::string zHeader = properties["z" + headerKey];
			zCol = findHeader(zHeader);
		}

		if (xCol < 0 || yCol < 0)
		{
			std::cout << "Cannot draw graph. Headers incompatible: " << xHeader << ", " << yHeader << std::endl;
			continue;
		}

		std::string minKey = "Min" + i_to_str(count);
		std::string maxKey = "Max" + i_to_str(count);
		std::string roundKey = "round" + i_to_str(count);
		double minX = -FLT_MAX;
		double minY = -FLT_MAX;
		double minZ = -FLT_MAX;
		double maxX = FLT_MAX;
		double maxY = FLT_MAX;
		double maxZ = FLT_MAX;
		bool round = (properties.count(roundKey) > 0);

		if (properties.count("x" + minKey)) minX = atof(properties["x" + minKey].c_str());
		if (properties.count("y" + minKey)) minY = atof(properties["y" + minKey].c_str());
		if (properties.count("z" + minKey)) minZ = atof(properties["z" + minKey].c_str());
		if (properties.count("x" + maxKey)) maxX = atof(properties["x" + maxKey].c_str());
		if (properties.count("y" + maxKey)) maxY = atof(properties["y" + maxKey].c_str());
		if (properties.count("z" + maxKey)) maxZ = atof(properties["z" + maxKey].c_str());
		
		if (minX == -FLT_MAX || maxX == FLT_MAX) minMaxCol(xCol, &minX, &maxX, false);
		if (minY == -FLT_MAX || maxY == FLT_MAX) minMaxCol(yCol, &minY, &maxY, round);

		if (hasZ)
		{
			if (minZ == -FLT_MAX || maxZ == FLT_MAX) 
			{
				minMaxCol(zCol, &minZ, &maxZ, round);
			}
		}

		_minZ = minZ;
		_maxZ = maxZ;

		if (count == 0)
		{
			// draw X axis
			double aValue = minX;
			double aStep = (maxX - minX) / (double)xIntervals;
			double xPrec = 3;
			if (maxX > 10) xPrec = 0;

			for (int i = 0; i < xIntervals + 1; i++)
			{
				double proportion = (double)i / (double)xIntervals;
				double pos = width * xAxis * proportion;
				png->drawLine(pos, 0, pos, 5, 0, 0, 0, 0);

				png->drawText(f_to_str(aValue, xPrec), pos, 30);

				aValue += aStep;
			}

			std::string xTitle = xHeader;
			if (properties.count("xTitle" + i_to_str(count)))
			{
				xTitle = properties["xTitle" + i_to_str(count)];
			}

			png->drawText(xTitle, width * xAxis * 0.5, 85);

			aValue = minY;
			aStep = (maxY - minY) / (double)yIntervals;
			double yPrec = 2;
			if (maxY > 10) yPrec = 0;

			for (int i = 0; i < yIntervals + 1; i++)
			{
				double proportion = (double)i / (double)yIntervals;
				double pos = -height * yAxis * proportion;
				png->drawLine(0, pos, -5, pos, 0, 0, 0, 0);
				png->drawText(f_to_str(aValue, yPrec), -50, pos);

				aValue += aStep;
			}

			std::string yTitle = yHeader;
			if (properties.count("yTitle" + i_to_str(count)))
			{
				yTitle = properties["yTitle" + i_to_str(count)];
			}

			png->drawText(yTitle, 70, -height * (yAxis + (1 - yAxis) / 4));
		}

		GraphStyle style = GraphStyleScatter;
		std::string styleKey = "style" + i_to_str(count);
		if (properties.count(styleKey))
		{
			if (properties[styleKey] == "scatter")
			{
				style = GraphStyleScatter;
			}
			else if (properties[styleKey] == "line")
			{
				style = GraphStyleLine;
			}
			else if (properties[styleKey] == "heatmap")
			{
				style = GraphStyleHeatMap;
			}
		}

		png_byte red = 0;
		png_byte blue = 0;
		png_byte green = 0;

		std::string colourKey = "colour" + i_to_str(count);
		if (properties.count(colourKey))
		{
			if (properties[colourKey] == "blue")
			{
				blue = 255;
			}
			if (properties[colourKey] == "red")
			{
				red = 255;
			}
			if (properties[colourKey] == "green")
			{
				green = 255;
			}
		}

		float transparency = 0.0;
		std::string transKey = "transparency" + i_to_str(count);
		if (properties.count(transKey))
		{
			transparency = atof(properties[transKey].c_str());
		}

		std::vector<vec2> points;

		for (int i = 0; i < entries.size(); i++)
		{
			if (entries[i][xCol] != entries[i][xCol] || entries[i][yCol] != entries[i][yCol])
			{
				continue;
			}

			points.push_back(make_vec2(entries[i][xCol], entries[i][yCol]));
		}

		if (style != GraphStyleScatter && style != GraphStyleHeatMap)
		{
			std::sort(points.begin(), points.end(), vec2_less_vec2);
		}

		if (points.size() == 0)
		{
			count++;
			continue;
		}

		double lastX = (points[0].x - minX) / (maxX - minX);
		double lastY = (points[0].y - minY) / (maxY - minY);
		double xBlockSize = xAxis * width / (double)fastStride;
		int yNumber = entries.size() / (double)fastStride;
		double yBlockSize = yAxis * height / yNumber;

		for (int i = 0; i < points.size(); i++)
		{
			double x = points[i].x;
			double y = points[i].y;

			double xProp = (x - minX) / (maxX - minX);
			double yProp = (y - minY) / (maxY - minY);

			if (xProp > 1 || yProp > 1 || xProp < 0 || yProp < 0)
			{
				continue;
			}

			xProp *= xAxis * width;
			yProp *= -yAxis * height;

			if (style == GraphStyleScatter)
			{
				const int length = 1;
				png->drawLine(xProp - length, yProp - length, xProp + length, yProp + length, transparency, red, green, blue);
				png->drawLine(xProp - length, yProp + length, xProp + length, yProp - length, transparency, red, green, blue);
			}
			else if (style == GraphStyleLine)
			{
				if (i > 0)
				{
					png->drawLine(lastX, lastY, xProp, yProp, transparency, red, green, blue);
				}
				lastX = xProp;
				lastY = yProp;
			}
			else if (style == GraphStyleHeatMap)
			{
				double zValue = entries[i][zCol];
				double propZ = (zValue - minZ) / (maxZ - minZ);
				
				if (propZ > 2) propZ = 2;

				if (propZ <= 0)
				{
					propZ = std::min(-propZ, 1.);
					red = 0;
					green = 0;
					blue = 255 - propZ * 255;
				}
				else if (propZ < 0.5)
				{
					/* we go blue. */
					propZ = (0.5 - propZ ) * 2.;
					red = 255 - propZ * 255;
					green = 255 - propZ * 255;
					blue = 255;
				}
				else if (propZ >= 1.0) /* We go red. */
				{
					propZ -= 1; 
					red = 255;
					green = propZ * 255;
					blue = 0;
				}
				else if (propZ >= 0.5) /* We go red. */
				{
					propZ = (propZ - 0.5) * 2.0;
					red = 255;
					green = 255 - propZ * 255;
					blue = 255 - propZ * 255;
				}

				for (int j = yProp - yBlockSize; j < yProp + 1; j++)
				{
					for (int k = xProp - xBlockSize; k < xProp + 1; k++)
					{
						png->setPixelColourRelative(k, j, red, green, blue);
					}
				}
			}
		}

		count++;
	}

	png->writeImageOutput();
}

std::string CSV::plotColumns(int col1, int col2)
{
	char *cols = std::getenv("COLUMNS");
	char *linesstr = getenv("LINES");

	//   std::cout << "Plotting " << entries.size() << " points" << std::endl;

	int columns = 100;

	if (cols && strlen(cols))
	{
		columns = atoi(cols);
	}

	int lines = 50;

	if (cols && strlen(linesstr))
	{
		lines = atoi(linesstr);
	}

	if (columns < 80) columns = 80;
	if (lines < 50) lines = 50;

	const int leftMargin = 8;
	const int bottomMargin = 4;

	int graphCols = columns - leftMargin;
	int graphLines = lines - bottomMargin;

	Plot plot;

	for (int i = 0; i < lines; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			plot[i][j] = PlotBlank;
		}
	}

	double xMin, xMax, yMin, yMax;

	if (!didSetMinMaxXY)
	{
		minMaxCol(col1, &xMin, &xMax);
		minMaxCol(col2, &yMin, &yMax);
	}
	else
	{
		xMin = minX;
		xMax = maxX;
		yMin = minY;
		yMax = maxY;
	}

	for (int i = 0; i < graphLines; i++)
	{
		plot[i][leftMargin] = PlotVerticalLine;
	}

	for (int i = 0; i < graphLines; i += ((double)graphLines / 10.))
	{
		plot[i][leftMargin - 1] = PlotHorizontalTickMark;
		std::ostringstream valueStream;
		double value = ((double)i / (double)graphLines) * (yMin - yMax) + yMax;

		valueStream << std::fixed << std::setprecision(2) << value;

		writeStringToPlot(valueStream.str(), &plot, i, 0);
	}

	for (int i = 0; i < graphCols; i++)
	{
		plot[graphLines][leftMargin + i] = PlotHorizontalLine;
	}

	writeStringToPlot(headers[col1], &plot, graphLines + bottomMargin - 1, graphCols / 2);

	for (int i = 0; i < graphCols; i += ((double)graphCols / 10.))
	{
		plot[graphLines + 1][leftMargin + i] = PlotVerticalTickMark;

		std::ostringstream valueStream;
		double value = ((double)i / (double)graphCols) * (xMax - xMin) + xMin;

		valueStream << std::fixed << std::setprecision(2) << value;

		writeStringToPlot(valueStream.str(), &plot, graphLines + 2, leftMargin + i - 2);
	}

	for (int i = 0; i < entries.size(); i++)
	{
		double xValue = entries[i][col1];
		double yValue = entries[i][col2];

		if (xValue < xMin || xValue > xMax || yValue < yMin || yValue > yMax)
		continue;

		double xProportion = (xValue - xMin) / (xMax - xMin);
		double yProportion = (yValue - yMin) / (yMax - yMin);

		int xCharsIn = xProportion * (double)(graphCols - 1);
		int yCharsIn = yProportion * (double)graphLines;

		int xFullCharsIn = xCharsIn + leftMargin + 2;
		int yFullCharsIn = (graphLines - yCharsIn);

		if (yFullCharsIn == graphLines)
		yFullCharsIn -= 1;

		char currentChar = plot[yFullCharsIn][xFullCharsIn];

		if (currentChar == ' ')
		{
			plot[yFullCharsIn][xFullCharsIn] = '.';
		}
		else if (currentChar == '.')
		{
			plot[yFullCharsIn][xFullCharsIn] = ':';
		}
		else if (currentChar == ':')
		{
			plot[yFullCharsIn][xFullCharsIn] = '^';
		}
		else if (currentChar == '^')
		{
			plot[yFullCharsIn][xFullCharsIn] = '*';
		}
		else if (currentChar == '*')
		{
			plot[yFullCharsIn][xFullCharsIn] = 'x';
		}
		else if (currentChar == 'x')
		{
			plot[yFullCharsIn][xFullCharsIn] = '#';
		}
		else if (currentChar == '#')
		{
			plot[yFullCharsIn][xFullCharsIn] = '%';
		}
		else if (currentChar == '%')
		{
			plot[yFullCharsIn][xFullCharsIn] = '@';
		}
	}

	std::ostringstream ascii;
	ascii << std::endl << headers[col2] << std::endl << std::endl;

	ascii << mapToAscii(plot);

	std::cout << ascii.str() << std::endl;

	return ascii.str();
}

void CSV::setValueForEntry(int entry, std::string header, double value)
{
	int column = findHeader(header);
	entries[entry][column] = value;
}

void CSV::resetColumn(std::string header, double value)
{
	int headerNum = findHeader(header);

	for (int i = 0; i < entries.size(); i++)
	{
		entries[i][headerNum] = value;
	}
}

CSV::~CSV()
{
	headers.clear();
	entries.clear();

	std::vector<std::string>().swap(headers);
	std::vector<Entry>().swap(entries);
}


