#include "codec.h"
#include "codec_handler.h"

//https://developers.google.com/streetview/publish/camm-spec

//TODO could tell from the working sample which types are used.
#include <set>
#include <iostream>
using namespace std;

namespace {
// Refactor step 39: CammHandler holds the full camm logic. CammHandler::
// search routes through the free gpmdSearch (preserving the pre-refactor
// quirk that camm's search dispatch went through gpmd, not through the
// long-dead Codec::cammSearch — which has now been removed).
struct CammHandler : CodecHandler {
	Match match(Codec &, const unsigned char *start, int maxlength) override {
		Match match;

		int lengths[] = { 12, 8, 12, 12, 12, 24, 14*4, 12 };
		while(start[0] == 0 && start[1] == 0 || start[3] == 0 && start[2] < 7) {
			int type = start[2];
			match.chances = 1<<20;
			int len = lengths[type] + 4;
			match.length += len;
			start += len;
			if(match.length > maxlength - 15*4)
				break;
		}
		return match;
	}
	Match search(Codec &, const unsigned char *start, int maxlength) override {
		return gpmdSearch(start, maxlength);
	}
};
} // namespace

std::unique_ptr<CodecHandler> makeCammHandler() {
	return std::unique_ptr<CodecHandler>(new CammHandler());
}
