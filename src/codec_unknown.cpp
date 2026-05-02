#include "codec.h"
#include "codec_handler.h"
#include "log.h"

////it might be encode each sample (1-4) bytes as a packet (AARGH!)

namespace {
// Refactor step 27: body inlined from the now-removed Codec::unknownMatch.
// Reads codec.name and codec.stats. Used as the final fallthrough handler
// for any codec name not matched by a more specific handler (and not pcm).
// Covers "rtmd", "samr", and any unrecognized codec.
struct UnknownHandler : CodecHandler {
	Match match(Codec &codec, const unsigned char *start, int maxlength) override {
		Match match;

		if(codec.name == "rtmd") {
			int32_t begin32 = readBE<int32_t>(start);

			if(codec.stats.beginnings32.count(begin32))
				match.chances = 10*codec.stats.beginnings32[begin32];
			match.length = 1024;
			return match;
		}
		//probably no samples
		if(!codec.stats.fixed_size)
			return match;

		if(codec.stats.fixed_size)
			match.length = codec.stats.fixed_size;

		if(codec.name == "samr") {
			//samr special thing might be 160... or different values.

			static const uint8_t amrnb_packed_size[16] = {
				13, 14, 16, 18, 20, 21, 27, 32, 6, 1, 1, 1, 1, 1, 1, 1
			};
			static const uint8_t amrwb_packed_size[16] = {
				18, 24, 33, 37, 41, 47, 51, 59, 61, 6, 1, 1, 1, 1, 1, 1
			};
			int mode = (start[0] >> 3) & 0xf;
			if (mode > 9 || (start[0] & 0x4) != 0x4) {
				match.chances = 0;
				return match;
			}
			match.length = amrnb_packed_size[mode];
			match.chances = 4;
		}




		int64_t begin64 = readBE<int64_t>(start);
		int32_t begin32 = readBE<int32_t>(start);

		/*if(stats.fixed_begin64) {
			if(stats.fixed_begin64 == begin64)
				match.chances = 1e40f;
			else
				match.chances = 0.0f;
			return match;
		}

		if(stats.fixed_begin32) {
			if(stats.fixed_begin32 == begin32)
				match.chances = 1e20f;
			else
				match.chances = 0.0f;
			return match;
		}*/

		//if we have a probability with the beginnig we can do something otherwise just lower it below 2 (for pcm)
		//we will have to look up for patterns, and when chances are really low, use search!
		if(codec.stats.beginnings32.count(begin32))
			match.chances = codec.stats.beginnings32[begin32];
		if(codec.stats.beginnings64.count(begin64))
			match.chances = codec.stats.beginnings64[begin64];

		//match.chances = 0.0f; //we really have no idea

		return match;
	}
};
} // namespace

std::unique_ptr<CodecHandler> makeUnknownHandler() {
	return std::unique_ptr<CodecHandler>(new UnknownHandler());
}
