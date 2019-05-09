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

#ifndef __vagabond__BaseParser_h__
#define __vagabond__BaseParser_h__

#include "BaseParser.h"

class Parser : public BaseParser
{
public:	
	Parser() : BaseParser() {};
	virtual ~Parser() {};

	ParserPtr shared_from_this()
	{
		return ToParserPtr(BaseParser::shared_from_this());
	}
	

	static ParserPtr processBlock(char *block);
};


#endif
