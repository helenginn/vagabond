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

#include "TextManager.h"
#include "font.h"
#include <cstdlib>

void TextManager::text_free(png_byte **pointer)
{
	free(*pointer);
}

void TextManager::text_malloc(png_byte **pointer, std::string text, int *width, int *height)
{
	if (!text.length())
	{
		return;
	}

	int squish = 3;
	int totalWidth = squish;
	int maxHeight = 0;

	for (int i = 0; i < text.length(); i++)
	{
		int charHeight = asciiDimensions[text[i]][0];
		int charWidth = asciiDimensions[text[i]][1];
		totalWidth += charWidth - squish;
		if (charHeight > maxHeight)
		{
			maxHeight = charHeight;
		}
	}

	totalWidth += squish;

	*width = totalWidth;
	*height = maxHeight;

	*pointer = (png_byte *)calloc(maxHeight * totalWidth, sizeof(png_byte));

	int currentX = 0;

	for (int i = 0; i < text.length(); i++)
	{
		char whichAscii = text[i];
		png_byte *chosenAscii = asciis[whichAscii];
		int chosenHeight = asciiDimensions[whichAscii][0];
		int chosenWidth = asciiDimensions[whichAscii][1];

		for (int k = 0; k < chosenHeight; k++)
		{
			for (int j = 0; j < chosenWidth; j++)
			{
				int xPos = currentX + j;
				int newPos = totalWidth * k + xPos;
				int asciiPos = chosenWidth * k + j;

				(*pointer)[newPos] += chosenAscii[asciiPos];
			}
		}

		currentX += chosenWidth - squish;
	}
}
