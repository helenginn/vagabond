// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#ifndef __vagabond__ParserTypes__
#define __vagabond__ParserTypes__

#include <string>
#include "vec3.h"
#include "mat3x3.h"

/** \cond SHOW_PARSER_PROPERTY_STRUCTS */

typedef void (*Encoder)(void *, void *, std::ostream &, int);
typedef char *(*Decoder)(void *, void *, char *block);

typedef struct
{
	std::string *stringPtr;
	std::string ptrName;
} StringProperty;

typedef struct
{
	double *doublePtr;
	std::string ptrName;
} DoubleProperty;

typedef struct
{
	mat3x3 *mat3x3Ptr;
	std::string ptrName;
} Mat3x3Property;

typedef struct
{
	vec3 *vec3Ptr;
	std::string ptrName;
} Vec3Property;

typedef struct
{
	std::vector<vec3> *vec3ArrayPtr;
	std::string ptrName;
} Vec3ArrayProperty;

typedef struct
{
	std::vector<mat3x3> *mat3x3ArrayPtr;
	std::string ptrName;
} Mat3x3ArrayProperty;

typedef struct
{
	bool *boolPtr; 
	std::string ptrName;
} BoolProperty;

typedef struct
{
	int *intPtr;
	std::string ptrName;
} IntProperty;

typedef struct
{
	void *objPtr; /** Object is passed to delegate for encoding/decoding */
	std::string ptrName;
	void *delegate; /** Decoder/Encoder of the delegate is called */
	Encoder encoder;
	Decoder decoder;
} CustomProperty;

/** \endcond */

#endif
