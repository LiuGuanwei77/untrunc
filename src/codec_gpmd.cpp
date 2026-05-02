#include "codec.h"
#include "codec_handler.h"

#include <set>
#include <iostream>
using namespace std;

namespace {
// Refactor step 29: gpmdMatch body inlined here. gpmdSearch is still
// kept as a Codec method for now because CammHandler::search also calls
// it (codec.gpmdSearch); a follow-up step will lift it to a free
// function so both handlers can share it.
struct GpmdHandler : CodecHandler {
	Match match(Codec &, const unsigned char *start, int maxlength) override {
		Match match;

		set<string> fourcc =  {"DEVC", "DVID", "DVNM", "STRM", "STNM", "RMRK",  "SCAL",
								 "SIUN", "UNIT", "TYPE", "TSMP", "TIMO", "EMPT"};
		if(!fourcc.count(string((char *)start, 4)))
			return match;

		match.chances = 1<<20;
		match.length = (readBE<int32_t>(start + 4) & 0xffff) + 8;
		return match;
	}
	Match search(Codec &codec, const unsigned char *start, int maxlength) override {
		return codec.gpmdSearch(start, maxlength);
	}
};
} // namespace

std::unique_ptr<CodecHandler> makeGpmdHandler() {
	return std::unique_ptr<CodecHandler>(new GpmdHandler());
}

Match Codec::gpmdSearch(const unsigned char *start, int maxlength) {
	Match match;
	const unsigned char *end = start + maxlength - 8;
	const  unsigned char *current = start;

	set<string> fourcc =  {"DEVC", "DVID", "DVNM", "STRM", "STNM", "RMRK",  "SCAL",
							 "SIUN", "UNIT", "TYPE", "TSMP", "TIMO", "EMPT"};

	while(current < end) {
		if(fourcc.count(string((char *)current, 4))) {
			match.chances = 1<<20;
			match.offset = current - start;
			return match;
		}
		current++;
	}
	return match;
}
