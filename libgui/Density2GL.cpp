//
//  Density2GL.cpp
//  VagabondViewer
//
//  Created by Helen Ginn on 22/7/2018.
//  Copyright © 2018 Ginn. All rights reserved.
//

#include "Density2GL.h"
#include "GLKeeper.h"
#include "../../libsrc/Crystal.h"
#include "../../libsrc/fftw3d.h"
#include <cstdarg>
#include <algorithm>

std::vector<GLuint> makeList(int count, ...)
{
	va_list arguments;
	va_start(arguments, count);

	std::vector<GLuint> list;

	for (int i = 0; i < count; i++)
	{
		GLuint val = va_arg(arguments, GLuint);
		list.push_back(val);
	}

	va_end(arguments);
	
	return list;
}

void initVertex(Vertex *vert)
{
	vert->pos[0] = 0;
	vert->pos[1] = 0;
	vert->pos[2] = 0;
	
	int random = rand() % 3;
	
	vert->color[0] = 0.8;
	vert->color[1] = 0.8;
	vert->color[2] = 0.8;
	vert->color[3] = 0.05;
}

void Density2GL::setupIndexTable()
{
	int y = _dims.x * 3;
	int z = _dims.y * _dims.x * 3;
	
	int y2 = y + 2;
	int x = 3;
	int x1 = 4;
	int x2 = 5;
	int xy2 = y + 5;
	int yz = y + z;
	int z1 = z + 1;
	int xz1 = z + 4;

	std::cout << y << ", ";
	_cubeIndices.resize(256);
	
	// type 1 checked
	_cubeIndices[0b00000001] = makeList(3, 0, 1, 2);
	_cubeIndices[0b11111110] = _cubeIndices[0b00000001];

	// type 1 checked
	_cubeIndices[0b00000010] = makeList(3, 0, x1, x2);
	_cubeIndices[0b11111101] = _cubeIndices[0b00000010];

	// type 2 checked
	_cubeIndices[0b00000011] = makeList(6, 1, x1, x2, 1, x2, 2);
	_cubeIndices[0b11111100] = _cubeIndices[0b00000011]; 

	// type 1 checked
	_cubeIndices[0b00000100] = makeList(3, 1, y, y2);
	_cubeIndices[0b11111011] = _cubeIndices[0b00000100];

	// type 2 checked
	_cubeIndices[0b00000101] = makeList(6, 0, y, y2, 0, y2, 2);
	_cubeIndices[0b11111010] = _cubeIndices[0b00000101];

	// type 3 checked
	_cubeIndices[0b00000110] = makeList(6, 0, x1, x2, 1, y, y2);
	_cubeIndices[0b11111001] = _cubeIndices[0b00000110];

	// type 4 checked
	_cubeIndices[0b00000111] = makeList(9, x2, 2, y2, x2, y2, x1,
	                                    x1, y2, y);
	_cubeIndices[0b11111000] = _cubeIndices[0b00000111];

	// type 1 checked
	_cubeIndices[0b00001000] = makeList(3, 2, z, z1);
	_cubeIndices[0b11110111] = _cubeIndices[0b00001000];

	// type 2 checked
	_cubeIndices[0b00001001] = makeList(6, 0, 1, z1, 0, z1, z);
	_cubeIndices[0b11110110] = _cubeIndices[0b00001001];

	// type 3 checked
	_cubeIndices[0b00001010] = makeList(6, 2, z, z1, 0, x1, x2);
	_cubeIndices[0b11110101] = _cubeIndices[0b00001010];

	// type 4 checked
	_cubeIndices[0b00001011] = makeList(9, x1, 1, z1, x1, z1, x2,
	                                    x2, x1, z);
	_cubeIndices[0b11110100] = _cubeIndices[0b00001011];

	// type 3 checked
	_cubeIndices[0b00001100] = makeList(6, 2, z, z1, 1, y, y2);
	_cubeIndices[0b11110011] = _cubeIndices[0b00001100];

	// type 4 fixed
	_cubeIndices[0b00001101] = makeList(9, 0, y, z, y, z, z1,
	                                    y, z1, y2);
	_cubeIndices[0b11110010] = _cubeIndices[0b00001101];

	// type 12 checked
	_cubeIndices[0b00001110] = makeList(9, 2, z, z1, 1, y, y2, 0, x1, x2);
	_cubeIndices[0b11110001] = _cubeIndices[0b00001110];

	// type 8 fixed again
	_cubeIndices[0b11110000] = makeList(12, x1, y, y2, x1, y2, x2,
	                                    y2, x2, z1, x2, z1, z);
	// type 8 fixed again
	_cubeIndices[0b00001111] = makeList(12, x1, y, y2, x1, y2, x2,
	                                    y2, x2, z1, x2, z1, z);

	// type 1 checked
	_cubeIndices[0b00010000] = makeList(3, x1, y, xy2);
	_cubeIndices[0b11101111] = _cubeIndices[0b00010000];

	// type 3 checked
	_cubeIndices[0b00010001] = makeList(6, 0, 1, 2, x1, y, xy2);
	_cubeIndices[0b11101110] = _cubeIndices[0b00010001];

	// type 2 checked
	_cubeIndices[0b00010010] = makeList(6, 0, y, x2, y, x2, xy2);
	_cubeIndices[0b11101101] = _cubeIndices[0b00010010];

	// type 4 checked
	_cubeIndices[0b00010011] = makeList(9, x2, 2, xy2, xy2, 2, 1, xy2, 1, y);
	_cubeIndices[0b11101100] = _cubeIndices[0b00010011];

	// type 2 checked
	_cubeIndices[0b00010100] = makeList(6, 1, x1, y2, x1, y2, xy2);
	_cubeIndices[0b11101011] = _cubeIndices[0b00010100];
	// 40

	// type 4 checked
	_cubeIndices[0b00010101] = makeList(9, 2, y2, xy2, xy2, 2, 0, 0, xy2, x1);
	_cubeIndices[0b11101010] = _cubeIndices[0b00010101];

	// type 4 checked
	_cubeIndices[0b00010110] = makeList(9, x2, xy2, y2, x, y2, 0, y2, 0, 1);
	_cubeIndices[0b11101001] = _cubeIndices[0b00010110];

	// type 5 checked again
	_cubeIndices[0b00010111] = makeList(6, 2, x2, xy2, xy2, 2, y2);
	// type 5 checked again
	_cubeIndices[0b11101000] = _cubeIndices[0b00010111];

	// type 10 checked
	_cubeIndices[0b00011000] = makeList(6, x1, y, xy2, 2, z, z1);
	_cubeIndices[0b11100111] = _cubeIndices[0b00011000];

	// type 11 fixed
	_cubeIndices[0b00011001] = makeList(9, 0, 1, z1, 0, z1, z, x1, y, xy2);
	_cubeIndices[0b11100110] = _cubeIndices[0b00011001];
	// 50

	// type 11 fixed
	_cubeIndices[0b00011010] = makeList(9, 2, z, z1, 0, y, x2, 
	                                    y, x2, xy2);
	_cubeIndices[0b11100101] = _cubeIndices[0b00011010];

	// type 9 checked
	_cubeIndices[0b00011011] = makeList(12, 1, y, xy2, 1, 2, z,
	                                    xy2, 1, z, xy2, z, x2);

	// type 9 checked
	_cubeIndices[0b11100100] = makeList(12, 1, y, xy2, xy2, xz1, z,
	                                    1, xy2, z, 1, z, z1);
	
	// type 11 checked
	_cubeIndices[0b00011100] = makeList(9, 1, x1, y2, x1, y2, xy2,
	                                    2, z, z1);
	                                    
	_cubeIndices[0b11100011] = _cubeIndices[0b00011100];

	// type 9 checked
	_cubeIndices[0b00011101] = makeList(12, x1, y2, xy2, x1, y2, z,
	                                    y2, z, z1, x1, z, 0);

	// type 9 checked
	_cubeIndices[0b11100010] = makeList(12, 0, z, x1, y2, z, z1,
	                                    xy2, y2, x1, y2, z, x1);
	
	// type 6 checked
	_cubeIndices[0b00011110] = makeList(12, x2, xy2, y2, x2, y2,
	                                     0, 0, y2, 1, 2, z, z1);
        
    // type 6 fixed
    _cubeIndices[0b11100001] = makeList(12, x2, xy2, y2, x2, y2, z,
                                        y2, z, z1, 0, 1, 2);

    // type 4 checked
	_cubeIndices[0b00011111] = makeList(9, x2, xy2, y2, x2, y2, z,
	                                    y2, z, z1);
	                                    
	_cubeIndices[0b11100000] = _cubeIndices[0b00011111];

	// type 1 checked
	_cubeIndices[0b00100000] = makeList(3, x2, xz1, z);
	_cubeIndices[0b11011111] = _cubeIndices[0b00100000];

	// type 1 checked
	_cubeIndices[0b00100001] = makeList(3, 0, 1, 2, x2, xz1, z);
	_cubeIndices[0b11011110] = _cubeIndices[0b00100001];

	// type 2 checked
	_cubeIndices[0b00100010] = makeList(6, 0, z, x1, x1, z, xz1);
	_cubeIndices[0b11011101] = _cubeIndices[0b00100010];

	// type 4 checked
	_cubeIndices[0b00100011] = makeList(9, 1, x1, xz1, xz1, 1, 2,
	                                    xz1, 2, z);
	_cubeIndices[0b11011100] = _cubeIndices[0b00100011];
	// 60

	// type 10 checked
	_cubeIndices[0b00100100] = makeList(6, x2, xz1, z, 1, y, y2);
	_cubeIndices[0b11011011] = _cubeIndices[0b00100100];

	// type 11 checked
	_cubeIndices[0b00100101] = makeList(9, x2, xz1, z, 0, y, y2,
	                                    0, y2, 2);
	_cubeIndices[0b11011010] = _cubeIndices[0b00100101];

	// type 11 checked
	_cubeIndices[0b00100110] = makeList(9, 1, y, y2, 0, z, x1,
	                                    x1, z, xz1);
	_cubeIndices[0b11011001] = _cubeIndices[0b00100110];

	// type 15 fixed
	_cubeIndices[0b00100111] = makeList(12, z, 2, y2, z, x1, y2,
	                                    z, x1, xz1, y2, 1, x1);

	// type 15 fixed
	_cubeIndices[0b11011000] = makeList(12, y, x1, y2, y2, x1, z,
	                                    z, x1, xz1, y2, z, 2); 

	// type 2 checked
	_cubeIndices[0b00101000] = makeList(6, x2, xz1, 2, 2, xz1, z1);

	_cubeIndices[0b11010111] = _cubeIndices[0b00101000];
	// 70

	// type 4 checked
	_cubeIndices[0b00101001] = makeList(9, 1, z1, xz1, xz1, 1, 0,
	                                    xz1, 0, x2); 

	_cubeIndices[0b11010110] = _cubeIndices[0b00101001];

	// type 4 checked
	_cubeIndices[0b00101010] = makeList(9, x1, xz1, z1, x1, z1, 0,
	                                    0, z1, 2);

	_cubeIndices[0b11010101] = _cubeIndices[0b00101010];

	// type 5 checked again
	_cubeIndices[0b00101011] = makeList(6, x1, xz1, z1, x1, z1, 1);

	// type 5 checked again
	_cubeIndices[0b11010100] = _cubeIndices[0b00101011];

	// type 11 checked
	_cubeIndices[0b00101100] = makeList(9, 1, y, y2, x2, xz1, 2, 2,
	                                    xz1, z1);
	_cubeIndices[0b11010011] = _cubeIndices[0b00101100];

	// type 9 checked
	_cubeIndices[0b00101101] = makeList(12, 0, xz1, x2, xz1, 0, y2,
	                                    0, 1, y2, xz1, y2, z1);

	// type 9 checked
	_cubeIndices[0b11010010] = makeList(12, 0, x2, xz1, 0, xz1, y2,
	                                    xz1, y2, yz, 0, y2, y);
	// 80

	// type 7 checked
	_cubeIndices[0b00101110] = makeList(12, 1, y, y2, x1, xz1, z1,
	                                    x1, z1, 0, 0, z1, 2);
	// type 7 checked
	_cubeIndices[0b11010001] = makeList(12, 0, 1, 2, x1, xz1, z1,
	                                    x1, z1, y, z1, y, y2);

	// type 4 checked
	_cubeIndices[0b00101111] = makeList(9, x1, xz1, z1, x1, z1, y,
	                                    z1, y, y2);
	_cubeIndices[0b11010000] = _cubeIndices[0b00101111];

	// type 3 checked
	_cubeIndices[0b00110000] = makeList(6, x1, y, xy2, x2, xz1, z);
	_cubeIndices[0b11001111] = _cubeIndices[0b00110000];

	// type 12
	_cubeIndices[0b00110001] = makeList(9, x1, y, xy2, x2, xz1, z,
	                                    0, 1, 2);
	_cubeIndices[0b11001110] = _cubeIndices[0b00110001];

	// type 4 checked
	_cubeIndices[0b00110010] = makeList(9, 0, y, z, z, y, xy2, y,
	                                    xy2, xz1);
	_cubeIndices[0b11001101] = _cubeIndices[0b00110010];
	// 90

	// type 8 checked again and again
	_cubeIndices[0b00110011] = makeList(12, 1, y, 2, y, 2, xy2, xy2,
	                                    2, z, xy2, z, xz1);

	// type 8 fixed again and again
	_cubeIndices[0b11001100] = makeList(12, y, 1, 2, 2, y, xy2, 2,
	                                    xy2, z, z, xy2, xz1); 
	
	// type 11
	_cubeIndices[0b00110100] = makeList(9, x2, xz1, z, 1, x1, y2,
	                                    x1, y2, xy2);
	_cubeIndices[0b11001011] = _cubeIndices[0b00110100];

	// type 6
	_cubeIndices[0b00110101] = makeList(9, x2, xz1, z, 2, y2, xy2,
	                                    xy2, 2, 0, 0, xy2, x1);
	_cubeIndices[0b11001010] = _cubeIndices[0b00110101];

	// type 9 checked
	_cubeIndices[0b00110110] = makeList(12, z, xz1, xy2, z, xy2, 1,
	                                    1, xy2, y, z, 1, 0);
	// type 9 checked
	_cubeIndices[0b11001001] = makeList(12, 0, 1, y2, 0, y2, z1,
	                                    y2, xz1, yz, 0, xz1, x2);

	// type 4 checked
	_cubeIndices[0b00110111] = makeList(9, 2, y2, xy2, 2, xy2, z,
	                                    z, xy2, xz1);
	_cubeIndices[0b11001000] = _cubeIndices[0b00110111];

	// type 11
	_cubeIndices[0b00111000] = makeList(9, x1, y, xy2, x2, xz1, 2,
	                                    2, xz1, z1);
	_cubeIndices[0b11000111] = _cubeIndices[0b00111000];

	// type 6
	_cubeIndices[0b00111001] = makeList(12, x1, y, xy2, 1, z1, xz1,
	                                    xz1, 1, 0, xz1, 0, x2);
	_cubeIndices[0b11000110] = makeList(12, 1, xz1, z1, 1, y, xz1,
	                                    xz1, y, xy2, 0, x1, x2);

	// type 15 checked
	_cubeIndices[0b00111010] = makeList(12, 0, 2, y, 2, y, xz1,
	                                    xz1, 2, z1, y, xz1, xy2);
	// type 15 checked
	_cubeIndices[0b11000101] = makeList(12, z1, 2, 0, z1, 0, xy2,
	                                    xy2, 0, y, z1, xy2, xz1);

	// type 4 checked
	_cubeIndices[0b00111011] = makeList(9, 1, z1, xz1, xz1, 1, y,
	                                    xz1, y, xy2);
	// type 4 checked
	_cubeIndices[0b11000100] = _cubeIndices[0b00111011];

	// type 14
	_cubeIndices[0b00111100] = makeList(12, 1, x1, y2, x1, y2, xy2,
	                                    x2, xz1, 2, 2, xz1, z1);
	// type 14
	_cubeIndices[0b11000011] = makeList(12, y2, xy2, xz1, y2, xz1,
	                                    z1, 1, x1, x2, 1, x2, 2);
	// type 11
	_cubeIndices[0b00111101] = makeList(9, 0, x1, x2, y2, xy2, xz1,
	                                    xz1, y2, z1);
	_cubeIndices[0b11000010] = _cubeIndices[0b00111101];

	// type 11
	_cubeIndices[0b00111110] = makeList(9, 0, 1, 2, y2, xy2, xz1,
	                                    xz1, y2, z1);
	_cubeIndices[0b11000001] = _cubeIndices[0b00111110];

	// type 2 checked
	_cubeIndices[0b00111111] = makeList(6, y2, xy2, xz1,
	                                    xz1, y2, z1);
	_cubeIndices[0b11000000] = _cubeIndices[0b00111111];

	// type 1 checked
	_cubeIndices[0b01000000] = makeList(3, y2, yz, z1);
	_cubeIndices[0b10111111] = _cubeIndices[0b01000000];

	// type 3
	_cubeIndices[0b01000001] = makeList(6, 0, 1, 2, y2, yz, z1);
	_cubeIndices[0b10111110] = _cubeIndices[0b01000001];

	// type 10
	_cubeIndices[0b01000010] = makeList(6, y2, yz, z1, x1, x2, 0);
	_cubeIndices[0b10111101] = _cubeIndices[0b01000010];

	// type 11
	_cubeIndices[0b01000011] = makeList(6, 1, x1, x2, 1, x2, 2,
	                                    y2, yz, z1);
	_cubeIndices[0b10111100] = _cubeIndices[0b01000011];

	// type 2 checked
	_cubeIndices[0b01000100] = makeList(6, y, 1, z1, y, z1, yz);
	_cubeIndices[0b10111011] = _cubeIndices[0b01000100];

	// type 4
	_cubeIndices[0b01000101] = makeList(9, 0, y, yz, 0, yz, 2,
	                                    yz, 2, z1);
	_cubeIndices[0b10111010] = _cubeIndices[0b01000101];
	
	// type 11
	_cubeIndices[0b01000110] = makeList(9, y, 1, z1, y, z1, yz,
	                                    0, x1, x2);
	_cubeIndices[0b10111001] = _cubeIndices[0b01000110];

	// type 15 (type 9?)
	_cubeIndices[0b01000111] = makeList(0);
	// type 15 (type 9?)
	_cubeIndices[0b10111000] = makeList(0);

	// type 2 fixed
	_cubeIndices[0b01001000] = makeList(6, 2, y2, yz, 2, yz, z);
	_cubeIndices[0b10110111] = _cubeIndices[0b01001000];

	// type 6 fixed
	_cubeIndices[0b01001001] = makeList(9, 0, z, yz, 0, yz, 1, 
	                                    yz, 1, y2);
	_cubeIndices[0b10110110] = _cubeIndices[0b01001001];

	// type 11
	_cubeIndices[0b01001010] = makeList(9, y2, 2, z, y2, z, yz,
	                                    0, x1, x2);
	_cubeIndices[0b10110101] = _cubeIndices[0b01001010];

	// type 15 checked
	_cubeIndices[0b01001011] = makeList(12, y2, yz, 1, 1, yz, x2,
	                                    yz, x2, z, 1, x2, x1);
	// type 15 checked
	_cubeIndices[0b10110100] = makeList(12, yz, y2, 1, yz, 1, x2,
	                                    x2, 1, x1, yz, x2, z);

	// type 6 checked
	_cubeIndices[0b01001100] = makeList(9, y, yz, z, y, z, 1, 
	                                    1, z, 2);
	_cubeIndices[0b10110011] = _cubeIndices[0b01001100];

	// type 5 checked
	_cubeIndices[0b01001101] = makeList(6, 0, y, yz, yz, 0, z);

	// type 5 checked
	_cubeIndices[0b10110010] = makeList(6, 0, yz, y, 0, yz, z);

	// type 6 checked
	_cubeIndices[0b01001110] = makeList(12, y, yz, z, y, z, 1,
	                                    1, z, 2, 0, x1, x2);

	// type 6 checked
	_cubeIndices[0b10110001] = makeList(12, y, z, yz, y, z, x1,
	                                    x1, z, x2, 0, 1, 2);

	// type 4 checked
	_cubeIndices[0b01001111] = makeList(9, y, yz, z, y, z, x1, 
	                                    z, x1, x2);
	_cubeIndices[0b10110000] = _cubeIndices[0b01001111];

	// type 3
	_cubeIndices[0b01010000] = makeList(6, y2, yz, z1, x1, y, xy2);
	_cubeIndices[0b10101111] = _cubeIndices[0b01010000];

	// type 12
	_cubeIndices[0b01010001] = makeList(9, 0, 1, 2, yz, z1, y2,
	                                    y, xy2, x1);
	_cubeIndices[0b10101110] = _cubeIndices[0b01010001];

	// type 11 checked
	_cubeIndices[0b01010010] = makeList(9, yz, z1, y2, 0, y, xy2,
	                                    0, xy2, x2);
	_cubeIndices[0b10101101] = _cubeIndices[0b01010010];

	// type 6 checked
	_cubeIndices[0b01010011] = makeList(12, 2, x2, xy2, xy2, 2, 1,
	                                    xy2, 1, y, y2, z1, yz);
	// type 6 fixed
	_cubeIndices[0b10101100] = makeList(12, 1, y, y2, 1, y, y2,
	                                    2, xy2, z1, xy2, z1, yz);

	// type 4 fixed
	_cubeIndices[0b01010100] = makeList(9, x1, 1, z1, x1, z1, xy2,
	                                    z1, xy2, yz);
	_cubeIndices[0b10101011] = _cubeIndices[0b01010100];

	// type 8 fixed 
	_cubeIndices[0b01010101] = makeList(12, 0, x1, 2, x1, 2, xy2, 2,
	                                    xy2, z1, xy2, z1, yz);
	// type 8 fixed
	_cubeIndices[0b10101010] = makeList(12, 0, x1, 2, x1, 2, xy2, 2,
	                                    xy2, z1, xy2, z1, yz);

	// type 15 fixed
	_cubeIndices[0b01010110] = makeList(12, 0, 1, x2, 1, x2, yz,
	                                    xy2, yz, x2, 1, yz, z1);
	// type 15 fixed
	_cubeIndices[0b10101001] = makeList(12, 0, x2, 1, 1, x2, yz,
	                                    z1, 1, yz, x2, xy2, yz);

	// type 6 fixed
	_cubeIndices[0b01010111] = makeList(9, xy2, x2, 2, 2, xy2, z1,
	                                    xy2, z1, yz);
	_cubeIndices[0b10101000] = _cubeIndices[0b01010111];

	// type 11 fixed
	_cubeIndices[0b01011000] = makeList(9, x1, y, xy2, y2, yz, z,
	                                    z, y2, 2);
	_cubeIndices[0b10100111] = _cubeIndices[0b01011000];

	// type 6 checked
	_cubeIndices[0b01011001] = makeList(12, x1, y, xy2, 0, z, yz,
	                                    0, yz, 1, yz, 1, y2);
	// type 6 fixed
	_cubeIndices[0b10100110] = makeList(12, 1, y, y2, 0, z, yz, 0,
	                                    yz, x1, yz, x1, xy2);

	// type 14 FIXME
//	_cubeIndices[0b01011010] = makeList(12, 

	// type 14 FIXME
//	_cubeIndices[0b10100101] = makeList(12, 

	// type 11
	_cubeIndices[0b01011011] = makeList(9, x1, xz1, 1, xz1, 1, y2,
	                                    y2, xz1, yz);
	_cubeIndices[0b10100100] = _cubeIndices[0b01011011];

	// type 9 checked
	_cubeIndices[0b01011100] = makeList(12, x1, xy2, yz, x1, yz, 2,
	                                    2, yz, z1, x1, 2, 1);
	// type 9 checked
	_cubeIndices[0b10100011] = makeList(12, yz, xy2, x1, yz, x1, 2,
	                                    x1, 2, 0, yz, 2, z);

	// type 4 checked
	_cubeIndices[0b01011101] = makeList(9, 0, z, yz, 0, yz, x1, 
	                                    x1, yz, xy2);
	_cubeIndices[0b10100010] = _cubeIndices[0b01011101];

	// type 11 checked
	_cubeIndices[0b01011110] = makeList(9, x2, xy2, yz, x2, yz, z, 
	                                    0, 1, 2);
	_cubeIndices[0b10100001] = _cubeIndices[0b01011110];

	// type 2 checked
	_cubeIndices[0b01011111] = makeList(6, x2, xy2, yz, x2, yz, z);
	_cubeIndices[0b10100000] = _cubeIndices[0b01011111];

	// type 3
	_cubeIndices[0b01100000] = makeList(6, x2, z, xz1, z1, y2, yz);
	_cubeIndices[0b10011111] = _cubeIndices[0b01100000];

	// type 12
	_cubeIndices[0b01100001] = makeList(9, x2, z, xz1, z1, y2, yz,
	                                    0, 1, 2);
	_cubeIndices[0b10011110] = _cubeIndices[0b01100001];

	// type 11 checked
	_cubeIndices[0b01100010] = makeList(9, 0, z, x1, x1, z, xz1, 
	                                    y2, z1, yz);
	_cubeIndices[0b10011101] = _cubeIndices[0b01100010];

	// type 6 fixed
	_cubeIndices[0b01100011] = makeList(12, 1, x1, xz1, 1, xz1, y,
	                                    xz1, y, z, y2, z1, yz);
	_cubeIndices[0b10011100] = _cubeIndices[0b01100011];

	// type 11 checked
	_cubeIndices[0b01100100] = makeList(9, z, x2, xz1, 1, y, z1, 
	                                    z1, y, yz);
	_cubeIndices[0b10011011] = _cubeIndices[0b01100100];

	// type 6
	_cubeIndices[0b01100101] = makeList(12, x2, z, xz1, yz, y, 0,
	                                    yz, 0, 2, 2, yz, z1);
	_cubeIndices[0b10011010] = _cubeIndices[0b01100101];

	// type 13
	_cubeIndices[0b10011001] = makeList(12, 0, 1, z1, 0, z1, z,
	                                    x1, xz1, y, y, yz, xz1);
	// type 13
	_cubeIndices[0b01100110] = makeList(12, 1, y, yz, 1, yz, z1,
	                                    0, z, xz1, xz1, 0, x1);
	                                    
	// type 11 fixed
	_cubeIndices[0b01100111] = makeList(9, 2, z, z1, x1, y, xz1,
	                                    y, xz1, yz);
	_cubeIndices[0b10011000] = _cubeIndices[0b01100111];

	// type 4 checked
	_cubeIndices[0b01101000] = makeList(9, x2, 2, y2, x2, y2, yz,
	                                    x2, yz, xz1);
	_cubeIndices[0b10010111] = _cubeIndices[0b01101000];

	// type 8 checked not sq
	_cubeIndices[0b01101001] = makeList(12, 0, 1, x2, x2, 1, xz1, 
	                                    1, xz1, y2, xz1, y2, yz);
	_cubeIndices[0b10010110] = _cubeIndices[0b01101001];

	// type 9 checked
	_cubeIndices[0b01101010] = makeList(12, x1, 0, 2, x1, 2, yz,
	                                    2, yz, y2, x1, yz, xz1);
	// type 9 checked
	_cubeIndices[0b10010101] = makeList(12, 2, 0, x1, 2, x1, yz,
	                                    x1, yz, xy2, 2, yz, y2);

	// type 4 checked
	_cubeIndices[0b01101011] = makeList(9, 1, x1, xz1, 1, xz1, y2, 
	                                    xz1, y2, yz);
	_cubeIndices[0b10010100] = _cubeIndices[0b01101011];

	// type 14 FIXME
	_cubeIndices[0b01101100] = makeList(0);
	// type 14 FIXME
	_cubeIndices[0b10010011] = makeList(0);

	// type 4 checked
	_cubeIndices[0b01101101] = makeList(9, 0, y, yz, 0, yz, x2, 
	                                    x2, yz, xz1);
	_cubeIndices[0b10010010] = _cubeIndices[0b01101101];

	// type 11 checked
	_cubeIndices[0b01101110] = makeList(9, 0, 1, 2, x1, y, yz, 
	                                    x1, yz, xz1);
	_cubeIndices[0b10010001] = _cubeIndices[0b01101110];

	// type 2 checked
	_cubeIndices[0b01101111] = makeList(6, x1, y, yz, x1, yz, xz1);
	_cubeIndices[0b10010000] = _cubeIndices[0b01101111];

	// type 12
	_cubeIndices[0b01110000] = makeList(9, z1, y2, yz, xy2, y, x1, 
	                                    z, x2, xz1);
	_cubeIndices[0b10001111] = _cubeIndices[0b01110000];

	// type 7
	_cubeIndices[0b01110001] = makeList(12, z1, y2, yz, xy2, y, x1,
	                                    z, x2, xz1, 0, 1, 2);
	_cubeIndices[0b10001110] = _cubeIndices[0b01110001];

	// type 6
	_cubeIndices[0b01110010] = makeList(12, z1, y2, yz, 0, y, z, 
	                                    y, z, xz1, y, xz1, xy2);
	_cubeIndices[0b10001101] = makeList(12, 0, y, z, y, z, y2, z,
	                                    y2, z1, xy2, yz, xz1);

	// type 7
	_cubeIndices[0b01110011] = makeList(9, 1, y, y2, 2, z, z1,
	                                    xy2, xz1, yz);
	_cubeIndices[0b10001100] = _cubeIndices[0b01110011];

	// type 6
	_cubeIndices[0b01110100] = makeList(12, x2, z, xz1, x1, y, z1, 
	                                    x1, z1, xy2, z1, xy2, yz);
	// type 6
	_cubeIndices[0b10001011] = makeList(0);

	// type 7
	_cubeIndices[0b01110101] = makeList(9, x2, 0, x1, 2, z, z1,
	                                    xz1, xy2, yz);
	_cubeIndices[0b10001010] = _cubeIndices[0b01110101];

	// type 11
	_cubeIndices[0b01110110] = makeList(9, xy2, yz, xz1, 0, 1, z1,
	                                    0, z1, z);
	_cubeIndices[0b10001001] = _cubeIndices[0b01110110];

	// type 2
	_cubeIndices[0b01110111] = makeList(6, xz1, xy2, yz, z1, 2, z);
	_cubeIndices[0b10001000] = _cubeIndices[0b01110111];

	// type 6
	_cubeIndices[0b01111000] = makeList(12, y2, z1, yz, xz1, x2, z,
	                                    x1, y, xy2, 2, z, z1);
	// type 6
	_cubeIndices[0b10000111] = makeList(12, xy2, xz1, yz, x2, 2, y2,
	                                    x2, y2, y, x2, y, x1);

	// type 12
	_cubeIndices[0b01111001] = makeList(9, 0, x1, x2, 1, y, y2, 
	                                    xy2, xz1, yz);
	_cubeIndices[0b10000110] = _cubeIndices[0b01111001];

	// type 11
	_cubeIndices[0b01111010] = makeList(9, xy2, yz, xz1, 0, y, y2, 
	                                    0, y2, 2);
	_cubeIndices[0b10000101] = _cubeIndices[0b01111010];

	// type 3
	_cubeIndices[0b01111011] = makeList(6, xy2, xz1, yz, 1, y, y2);
	_cubeIndices[0b10000100] = _cubeIndices[0b01111011];

	// type 11
	_cubeIndices[0b01111100] = makeList(9, xy2, xz1, yz, 1, x1, x2, 
	                                    1, x2, 2);
	_cubeIndices[0b10000011] = _cubeIndices[0b01111100];

	// type 3
	_cubeIndices[0b01111101] = makeList(6, xy2, xz1, yz, 0, x1, x2);
	_cubeIndices[0b10000010] = _cubeIndices[0b01111101];
	
	// type 10
	_cubeIndices[0b01111110] = makeList(6, xy2, xz1, yz, 0, 1, 2);
	_cubeIndices[0b10000001] = _cubeIndices[0b01111110];

	// type 1 checked
	_cubeIndices[0b01111111] = makeList(3, xy2, xz1, yz);
	_cubeIndices[0b10000000] = _cubeIndices[0b01111111];
}

vec3 Density2GL::getCentreOffset()
{
	vec3 halfsize = make_vec3(_dims.x, _dims.y, _dims.z);
	vec3_mult(&halfsize, -_resolution * 0.5);
	vec3 central = vec3_add_vec3(_offset, halfsize);
	
	return central;
}

void Density2GL::makeUniformGrid()
{
	/* Find sensible dimensions eventually */
	_dims.x = 20;
	_dims.y = 20;
	_dims.z = 20;
	
	/* Make a series of vertices which are three times more
	 * numerous than the number of voxels. */
	
	int total = _dims.x * _dims.y * _dims.z;
	total *= 3;
	_vertices.resize(total);
	long c = 0;
	vec3 central = getCentreOffset();
	
	for (int z = 0; z < _dims.z; z++)
	{
		for (int y = 0; y < _dims.y; y++)
		{
			for (int x = 0; x < _dims.x; x++)
			{
				for (int i = 0; i < 3; i++)
				{
					vec3 xyz = make_vec3(x, y, z);
					vec3_mult(&xyz, _resolution);
					vec3_add_to_vec3(&xyz, central);

					double add0 = (i == 0) ? _resolution / 2 : 0;
					double add1 = (i == 1) ? _resolution / 2 : 0;
					double add2 = (i == 2) ? _resolution / 2 : 0;

					initVertex(&_vertices[c]);
					_vertices[c].pos[0] = xyz.x + add0;
					_vertices[c].pos[1] = xyz.y + add1;
					_vertices[c].pos[2] = xyz.z + add2;
					c++;
				}
			}
		}
	}
}

void Density2GL::calculateContouring(CrystalPtr crystal)
{
	const int shiftCount = 8;
	vec3 shifts[] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
	                 {1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};

	FFTPtr fft = crystal->getFFT();
	mat3x3 real2frac = crystal->getReal2Frac();
	
	long count = -3;
	size_t unhandled = 0;
	size_t handled = 0;
	vec3 central = getCentreOffset();

	for (int z = 0; z < _dims.z; z++)
	{
		for (int y = 0; y < _dims.y; y++)
		{
			for (int x = 0; x < _dims.x; x++)
			{
				int bit = 0;
				
				count += 3;
				
				if (x >= _dims.x - 1 || y >= _dims.y - 1 ||
				    z >= _dims.z - 1)
				{
					continue;
				}

				for (int i = 0; i < shiftCount; i++)
				{
					vec3 xyz = make_vec3(x, y, z);
					vec3_add_to_vec3(&xyz, shifts[i]);
					vec3_mult(&xyz, _resolution);
					vec3_add_to_vec3(&xyz, central);

					mat3x3_mult_vec(real2frac, &xyz);

					double value = fft->getRealFromFrac(xyz);
					int crosses = (value > _threshold);

					bit |= (crosses << i);
				}
				
				/* Bit does not require any mesh if fully over
  				 * or under threshold */
				if (bit == 0 || bit == 255)
				{
					continue;
				}
				
				std::vector<GLuint> someIndices = _cubeIndices[bit];
				
				if (someIndices.size() == 0) 
				{
					unhandled++;
					continue;
				}

				handled++;

				for (int i = 0; i < someIndices.size(); i++)
				{
					someIndices[i] += count;
				}
				
				_indices.reserve(_indices.size() + someIndices.size());
				_indices.insert(_indices.end(), someIndices.begin(),
				                someIndices.end());
			}
		}
	}
	
//	std::cout << "Handled " << handled << " out of " << handled + unhandled;
//	std::cout << " (" << (double)handled / (double)(handled + unhandled) * 100 << ")";
//	std::cout << std::endl;
}

void Density2GL::nudgeDensity(int dir)
{
	if (dir > 0)
	{
		_threshold += 0.1;
	}
	else if (dir < 0)
	{
		_threshold -= 0.1;
	}
	
	recalculate();
	makeNewDensity();
}

void Density2GL::render()
{
	if (!_visible)
	{
		return;
	}
	
	bool locked = _renderLock.try_lock();
	
	if (locked)
	{
		if (_recalculate)
		{
			rebindProgram();
		}

		reorderIndices();
		GLObject::render();
		_renderLock.unlock();
	}
	
	mat4x4 inv = mat4x4_inverse(modelMat);
	vec3 centre = _keeper->getCentre();
	vec3 pan = _keeper->getTranslation();
	pan.z = 0;
	vec3_mult(&pan, -1);
	vec3_add_to_vec3(&centre, pan);
	centre.x = 0;
	centre.y = 0;
	vec3 newPos = mat4x4_mult_vec(inv, centre);

	vec3 movement = vec3_subtract_vec3(newPos, _offset);
	
	if (vec3_length(movement) > 1)
	{
		_offset = newPos;
		_recalculate = true;
		makeNewDensity();
	}
}

void Density2GL::makeNewDensity(CrystalPtr crystal)
{
	if (!crystal && !_crystal)
	{
		return;
	}
	else if (!crystal)
	{
		crystal = _crystal;
	}
	
	_crystal = crystal;
	
	_renderLock.lock();
	clearVertices();
	
	makeUniformGrid();

	if (!_cubeIndices.size())
	{
		setupIndexTable();
	}
	calculateContouring(crystal);
	reorderIndices();
	recalculate();
	_renderLock.unlock();
}
