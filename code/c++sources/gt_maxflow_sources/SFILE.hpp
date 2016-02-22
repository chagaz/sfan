/**********************************************************************************************//**
 *
 *	\file		SFILE.hpp
 *	\brief		SFILE class definition.
 *	\author		Killian Poulaud
 *	\version	0.1
 *	\date		2015/08/07
 *
 *	SFILE is a class that allow to do operations in a c_string like in I/O files.
 *
 *************************************************************************************************/
 
#pragma once

#include <cstdlib>
#include <cstddef>

#ifndef nullptr
#define nullptr NULL
#endif



namespace KCS
{
	
	class SFILE
	{
	  public:
		typedef unsigned long long ptr_size_t;
		
	  private:
		char *ptr;
		ptr_size_t cnt;
		ptr_size_t pos;
		bool free_ptr;
		
	  public:
		SFILE(const char *ptr_ = nullptr, bool free_ptr_ = false);
		SFILE(const SFILE &obj);
		~SFILE();
		
		SFILE& open(const char *c_str);
		void close();
		
		int getc();
		char* gets(char *tr, int num);
		int scans(const char *format, ...);
		size_t read(void *ptr, size_t size, size_t count);
		
		ptr_size_t tell();
		int seek(ptr_size_t offset, int origin);
		int eos() const;
		
		const char* getPtrClone() const;
		ptr_size_t getCnt() const;
		ptr_size_t getPos() const;
		
		bool operator!() const;
		SFILE& operator=(const SFILE &obj);
	};
	
}