#ifndef __Pic2Header_h__
#define __Pic2Header_h__

#include <stdio.h>

typedef struct
{
	size_t height;
	size_t width;
	size_t bytes_per_pix;
	void *data;
} Picture;

#endif
