#include "codec.h"
#include "log.h"
#include <string.h>
//#include "avlog.h"

extern "C" {
#include <libavformat/avformat.h>
#include <libavcodec/avcodec.h>
}

#include <iostream>
using namespace std;


inline int get_utf8_length(const char* text, const char* text_end)
{
	const char* start = text;
	uint32_t c;
	int err = 0;
	c = text < text_end ? (uint8_t)*text++ : (err = 1, 0);
	{
		uint32_t top = (c & 128) >> 1;
		if ((c & 0xc0) == 0x80 || c >= 0xFE)
		{
			goto error;
		}
		while (c & top) {
			unsigned int tmp = (text < text_end ? (uint8_t)*text++ : (err = 1, 0)) - 128;
			if (tmp >> 6)
			{
				goto error;
			}
			c = (c << 6) + tmp;
			top <<= 5;
		}
		c &= (top << 1) - 1;
	}
	if (err)
		goto error;
	return (int)(text - start);
error:
	return -1;
}

Match Codec::subtitleMatch(const unsigned char* start, int maxlength) {

	if (!context)
		throw string("Missing context for tx3g codec.");

	Match match;

	//TODO find more sensible values for these maxes
	uint32_t max_text_length = std::min(maxlength, 4096);
	uint32_t max_atom_length = 4096; //?

	int length = 0;
	
	uint16_t size;
	readBE<uint16_t>(size, start);

	if (size == 0 && incomplete_text_chunk) {
		match.chances = 0.1;
		match.length = 2;
		return match;
	}
	//probably there is a max size here!
	if (size > max_text_length || size <= 2 || size > maxlength - 2)
		return match;

	length += 2;
	const char *text_start = reinterpret_cast<const char*>(start + 2);
	const char* text = text_start;
	const char* text_end = text + size;
	while (text < text_end)
	{
		int len = get_utf8_length(text, text_end);
		if (len < 0)
			return match;
		text += len;
	}
	length += (int)(text - text_start);
		
	match.chances = 10000;
	match.length = length;
	incomplete_text_chunk = true;
	
	return match;
}
