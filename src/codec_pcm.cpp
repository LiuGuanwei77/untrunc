#include "codec.h"
#include "codec_handler.h"
#include "log.h"

////it might be encode each sample (1-4) bytes as a packet (AARGH!)

namespace {
// Refactor step 26: body inlined from the now-removed Codec::pcmMatch.
// Reads codec.stats.fixed_size (the Codec member that the original
// Codec::pcmMatch reached through implicit `this`).
struct PcmHandler : CodecHandler {
	Match match(Codec &codec, const unsigned char *start, int maxlength) override {
		Match match;
		if(codec.stats.fixed_size)
			match.length = codec.stats.fixed_size;
		match.chances = 2.0f; //we really have no idea
		return match;
	}
};
} // namespace

std::unique_ptr<CodecHandler> makePcmHandler() {
	return std::unique_ptr<CodecHandler>(new PcmHandler());
}
