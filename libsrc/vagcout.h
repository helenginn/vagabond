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

#ifndef __vagabond__vagcout__
#define __vagabond__vagcout__

#include <fstream>
#include <iostream>
#include <streambuf>
#include "FileReader.h"
#include "Notifiable.h"
#include "shared_ptrs.h"

template <typename charT, typename traits = std::char_traits<charT> >
class vagcout : public std::basic_streambuf<charT, traits>
{
public:
	static const size_t BUFF_SIZE = 1024;
	typedef traits traits_type;
	typedef typename traits_type::int_type int_type;
	typedef typename traits_type::pos_type pos_type;
	typedef typename traits_type::off_type off_type;
	
	explicit vagcout(std::streambuf *buf) : out_buf_(new charT[BUFF_SIZE])
	{
		this->_notify = NULL;
		this->original_cout = buf;
		outfile.open(FileReader::addOutputDirectory("vagabond.log"));
		this->setp(out_buf_, out_buf_ + BUFF_SIZE - 1);
	}
	
	~vagcout()
	{
		delete[] out_buf_;
		std::cout.rdbuf(original_cout);
		outfile.flush();
		outfile.close();
	}
	
	void setNotify(Notifiable *notify)
	{
		_notify = notify;
	}
protected:
	virtual int_type overflow(int_type c)
	{
		charT *ibegin = this->out_buf_;
		charT *iend = this->pptr();

		this->setp(out_buf_, out_buf_ + BUFF_SIZE + 1);

		if (!traits_type::eq_int_type(c, traits_type::eof()))
		{
			*iend++ = traits_type::to_char_type(c); // ew
		}

		int_type ilen = iend - ibegin;

		std::cout.rdbuf(original_cout);
		out_buf_[ilen] = '\0';
		std::cout << out_buf_;
		outfile << out_buf_;
		outfile << std::flush;
		
		if (_notify)
		{
			_notify->appendToLog(out_buf_);
		}

		std::cout.rdbuf(this);

		return traits_type::not_eof(c);
	}


	virtual int_type sync()
	{
		return traits_type::eq_int_type(this->overflow(traits_type::eof()),
		                                traits_type::eof()) ? -1 : 0;
	}

private:
	charT *out_buf_;
	
	std::streambuf *original_cout;
	std::ofstream outfile;

	Notifiable *_notify;
};


#endif
