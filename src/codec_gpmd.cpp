#include "codec.h"
#include "codec_handler.h"

#include <set>
#include <iostream>
using namespace std;

// Refactor step 30: lifted out of Codec::gpmdSearch (now removed). Body
// is byte-identical. Shared by GpmdHandler::search and CammHandler::search.
Match gpmdSearch(const unsigned char *start, int maxlength) {
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

namespace {
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
	Match search(Codec &, const unsigned char *start, int maxlength) override {
		return gpmdSearch(start, maxlength);
	}
};
} // namespace

std::unique_ptr<CodecHandler> makeGpmdHandler() {
	return std::unique_ptr<CodecHandler>(new GpmdHandler());
}
