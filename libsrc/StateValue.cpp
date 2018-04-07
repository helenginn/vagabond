// vagabond
//
// Created by Helen Ginn
// Copyright (c) 2018 Helen Ginn. All rights reserved.

#include "StateValue.h"
#include "Shouter.h"
#include "Parser.h"
#include "Bond.h"

void StateValue::applyToParser(Parser *parser)
{
	bool found = false;

	for (int i = 0; i < parser->_stringProperties.size(); i++)
	{
		if (parser->_stringProperties[i].ptrName == _ptrName)
		{
			*parser->_stringProperties[i].stringPtr = _string;

			found = true;
		}
	}

	for (int i = 0; i < parser->_doubleProperties.size(); i++)
	{
		if (parser->_doubleProperties[i].ptrName == _ptrName)
		{
			*parser->_doubleProperties[i].doublePtr = _double;

			found = true;
		}
	}

	for (int i = 0; i < parser->_mat3x3Properties.size(); i++)
	{
		if (parser->_mat3x3Properties[i].ptrName == _ptrName)
		{
			*parser->_mat3x3Properties[i].mat3x3Ptr = _mat3x3;

			found = true;
		}
	}

	for (int i = 0; i < parser->_vec3Properties.size(); i++)
	{
		if (parser->_vec3Properties[i].ptrName == _ptrName)
		{
			*parser->_vec3Properties[i].vec3Ptr = _vec3;

			found = true;
		}
	}

	for (int i = 0; i < parser->_vec3ArrayProperties.size(); i++)
	{
		if (parser->_vec3ArrayProperties[i].ptrName == _ptrName)
		{
			*parser->_vec3ArrayProperties[i].vec3ArrayPtr = _vec3Array;

			found = true;
		}
	}

	for (int i = 0; i < parser->_mat3x3ArrayProperties.size(); i++)
	{
		if (parser->_mat3x3ArrayProperties[i].ptrName == _ptrName)
		{
			*parser->_mat3x3ArrayProperties[i].mat3x3ArrayPtr = _mat3x3Array;

			found = true;
		}
	}

	for (int i = 0; i < parser->_boolProperties.size(); i++)
	{
		if (parser->_boolProperties[i].ptrName == _ptrName)
		{
			*parser->_boolProperties[i].boolPtr = _bool;

			found = true;
		}
	}

	for (int i = 0; i < parser->_intProperties.size(); i++)
	{
		if (parser->_intProperties[i].ptrName == _ptrName)
		{
			*parser->_intProperties[i].intPtr = _int;

			found = true;
		}
	}

	for (int i = 0; i < parser->_customProperties.size(); i++)
	{
		if (parser->_customProperties[i].ptrName == _ptrName)
		{
			CustomProperty property = parser->_customProperties[i];
			Decoder decoder = property.decoder;
			void *delegate = property.delegate;
			void *ptr = property.objPtr;
			char *custom = (char *)_custom.c_str();
			incrementIndent(&custom);
			(*decoder)(delegate, ptr, custom);
			repairCustom();
			
			static_cast<Bond *>(delegate)->postParseTidy();

			found = true;
		}
	}

	if (!found)
	{
		shout_at_helen("Missing property of type " + _ptrName
		               + " for parser " + parser->getAbsolutePath());	
	}
}

void StateValue::repairCustom()
{
	for (int i = 0; i < _custom.length(); i++)
	{
		if (_custom[i] == '\0')
		{
			_custom[i] = ' ';
		}	
	}	
}
