#include "codec.h"
#include "codec_handler.h"

#include <set>
#include <iostream>
using namespace std;

namespace {
// Refactor step 9: thin handler that delegates to Codec::gpmdMatch and
// Codec::gpmdSearch. The "camm" codec also dispatches search through
// gpmdSearch — that quirk will be preserved when camm is migrated.
struct GpmdHandler : CodecHandler {
	Match match(Codec &codec, const unsigned char *start, int maxlength) override {
		return codec.gpmdMatch(start, maxlength);
	}
	Match search(Codec &codec, const unsigned char *start, int maxlength) override {
		return codec.gpmdSearch(start, maxlength);
	}
};
} // namespace

std::unique_ptr<CodecHandler> makeGpmdHandler() {
	return std::unique_ptr<CodecHandler>(new GpmdHandler());
}
Match Codec::gpmdMatch(const unsigned char *start, int maxlength) {
	Match match;

	set<string> fourcc =  {"DEVC", "DVID", "DVNM", "STRM", "STNM", "RMRK",  "SCAL",
							 "SIUN", "UNIT", "TYPE", "TSMP", "TIMO", "EMPT"};
	if(!fourcc.count(string((char *)start, 4)))
		return match;

	match.chances = 1<<20;
	match.length = (readBE<int32_t>(start + 4) & 0xffff) + 8;
	return match;
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

