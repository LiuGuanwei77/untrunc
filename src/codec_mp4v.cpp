#include "codec.h"
#include "log.h"
//#include "avlog.h"
extern "C" {
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
}

using namespace std;

Match Codec::mp4vSearch(const unsigned char *start, int maxlength) {
	Match match;
	for(int offset = 0; offset < maxlength - 8; offset++) {
		int32_t begin32 = readBE<int32_t>(start + offset);
		if(begin32 == 0x1b3 || begin32 == 0x1b6) {
			match.offset = offset;
			match.chances = 1<<20;
			break;
		}
	}
	return match;
}

Match Codec::mp4vMatch(const unsigned char *start, int maxlength) {

	Match match;
	if(!context)
		return match;

	int32_t begin32 = readBE<int32_t>(start);
	if(begin32 != 0x1b3 && begin32 != 0x1b6)
		return match;
	match.chances = 1<<20;


	uint32_t duration = 0;

	int consumed = -1;

	{
		//AvLog useAvLog();
		av_log_set_level(0);

		static AVPacket* packet = av_packet_alloc();
		static AVFrame* frame = av_frame_alloc();

		packet->data = const_cast<unsigned char*>(start);
		packet->size = maxlength;
		int got_frame = 0;

		consumed = context->codec->decode(context, frame, &got_frame, packet);

//		bool keyframe = frame->key_frame;
//		not a frame? = !got_frame;
	}

	if(consumed == maxlength) {
		Log::debug << "Codec can't determine length of the packet.";
		match.chances = 4;
		consumed = 0; //unknown length
	} else if(consumed < 0) {
		match.chances = 0.0f;
		consumed = 0;
	}

	match.duration = duration;
	match.length = consumed;
	return match;
}
