#include "codec.h"
#include "codec_handler.h"

//https://developers.google.com/streetview/publish/camm-spec

//TODO could tell from the working sample which types are used.
#include <set>
#include <iostream>
using namespace std;

namespace {
// Refactor step 10: thin handler that delegates to Codec::cammMatch and,
// for search, Codec::gpmdSearch. The pre-refactor Codec::search dispatcher
// also routed camm's search to gpmdSearch (Codec::cammSearch exists but
// was never dispatched), so we preserve that quirk here.
struct CammHandler : CodecHandler {
	Match match(Codec &codec, const unsigned char *start, int maxlength) override {
		return codec.cammMatch(start, maxlength);
	}
	Match search(Codec &codec, const unsigned char *start, int maxlength) override {
		return codec.gpmdSearch(start, maxlength);
	}
};
} // namespace

std::unique_ptr<CodecHandler> makeCammHandler() {
	return std::unique_ptr<CodecHandler>(new CammHandler());
}
Match Codec::cammMatch(const unsigned char *start, int maxlength) {
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

Match Codec::cammSearch(const unsigned char *start, int maxlength) {
	Match match;
	const unsigned char *end = start + maxlength - 8;
	const  unsigned char *current = start;


	while(current < end) {
		if(start[0] == 0 && start[1] == 0 || start[3] == 0 && start[2] < 7) {
			match.chances = 1<<20;
			match.offset = current - start;
			return match;
		}
		current++;
	}
	return match;
}

