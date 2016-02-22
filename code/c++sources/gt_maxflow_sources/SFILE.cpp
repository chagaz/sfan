/**********************************************************************************************//**
 *
 *	\file		SFILE.cpp
 *	\brief		SFILE class methodes definition.
 *	\author		Killian Poulaud
 *	\version	0.1
 *	\date		2015/08/07
 *
 *	SFILE is a class that allow to do operations in a c_string like in I/O files.
 *
 *************************************************************************************************/
 
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cstddef>
#include <cstdarg>
#include <limits>
#include <algorithm>
#include <new>

#include "SFILE.hpp"

using namespace std;
using namespace KCS;



SFILE::SFILE(const char *ptr_, bool free_ptr_) :
	ptr(nullptr),
	cnt(0),
	pos(0),
	free_ptr(free_ptr_)
{
	open(ptr_);
}



SFILE::SFILE(const SFILE &obj) :
	SFILE(obj.getPtrClone(), true)
{
}



SFILE::~SFILE()
{
	close();
}



SFILE& SFILE::open(const char *ptr_)
{
	ptr = const_cast<char*>(ptr_);
	
	for (cnt = 0; cnt < numeric_limits<ptr_size_t>::max(); cnt++)
	{
		if (ptr[cnt] == '\0')
			break;
	}
	
	return *this;
}



void SFILE::close()
{
	if (free_ptr == true && ptr != nullptr)
		delete []ptr;
	
	ptr = nullptr;
	cnt = 0;
	pos = 0;
	free_ptr = false;
}



int SFILE::getc()
{
	if (!eos())
	{
		pos++;
		return ptr[pos - 1];
	}
	
	else
		return EOF;
}



char* SFILE::gets(char *tr, int num)
{
	int i;
	for (i = 0; (i < num - 1) & !eos(); i++)
	{	
		tr[i] = ptr[pos];
		pos++;
		
		if (tr[i] == 10)
			break;
	}
	
	if (i > 0)
	{
		tr[i] = '\0';
		return tr;
	}
	
	else
		return nullptr;
}



int SFILE::scans(const char *format, ...)
{
	va_list arg;
	int done;
	short token = 0;

	va_start(arg, format);
	done = vsscanf((ptr + pos), format, arg);
	va_end(arg);
	
	if (done > 0)
	{	
		for (int match_count = 0; (match_count < done) & !eos(); pos++)
		{
			if (ptr[pos] > 32)
				token = 1;
				
			else
			{
				if (token == 1)
					match_count++;
				
				token = 0;
			}
		}
	}

	return done;
}



size_t SFILE::read(void *buf, size_t size, size_t count)
{
	size_t bytes_requested = size * count;
	size_t bytes_read = 0;
	
	if (bytes_requested == 0)
		return 0;
	
	for (size_t i = 0; (bytes_read < bytes_requested) & !eos(); i++)
	{
		((char*)buf)[i] = ptr[pos];
		pos++;
		bytes_read++;
	}
	
	return (bytes_requested == bytes_read) ? count : bytes_read / size;
}



SFILE::ptr_size_t SFILE::tell()
{
	return pos;
}



int SFILE::seek(ptr_size_t offset, int origin)
{
	switch(origin)
	{
		case SEEK_SET:
			pos = 0;
			break;
			
		case SEEK_END:
			pos = cnt;
			break;
			
		case SEEK_CUR:
			;
	}
	
	pos += offset;
	
	if (pos > cnt)
		pos = cnt;
	
	return 0;
}



int SFILE::eos() const
{
	return (pos < cnt) ? 0 : 1;
}



const char* SFILE::getPtrClone() const
{
	char *new_ptr;
	
	try
	{
		int ptr_length = strlen(ptr);
		
		new_ptr = new char [ptr_length + 1];
		fill(new_ptr, new_ptr + ptr_length, '\0');
		
		strncpy(new_ptr, ptr, ptr_length);
		
		return new_ptr;
	}
	catch (bad_alloc &ba)
	{
		return nullptr;
	}
}



SFILE::ptr_size_t SFILE::getCnt() const
{
	return cnt;
}



SFILE::ptr_size_t SFILE::getPos() const
{
	return pos;
}



bool SFILE::operator!() const
{
	return (ptr == nullptr) ? true : false;
}



SFILE& SFILE::operator=(const SFILE &obj)
{
	ptr = const_cast<char*>(obj.getPtrClone());
	cnt = obj.getCnt();
	pos = obj.getPos();
	free_ptr = true;
	
	return *this;
}