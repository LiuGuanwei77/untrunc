#include "codec.h"
#include "codec_handler.h"
#include "log.h"

////it might be encode each sample (1-4) bytes as a packet (AARGH!)

namespace {
// Refactor step 18: thin handler that delegates to Codec::pcmMatch.
// Installed in Codec::parse() when the pcm flag is set (i.e. for codec
// names raw/twos/sowt/in24/in32/fl32/fl64/alaw/lpcm/ulaw).
struct PcmHandler : CodecHandler {
	Match match(Codec &codec, const unsigned char *start, int maxlength) override {
		return codec.pcmMatch(start, maxlength);
	}
};
} // namespace

std::unique_ptr<CodecHandler> makePcmHandler() {
	return std::unique_ptr<CodecHandler>(new PcmHandler());
}

Match Codec::pcmMatch(const unsigned char *start, int maxlength) {

	Match match;
	if(stats.fixed_size)
		match.length = stats.fixed_size;
	match.chances = 2.0f; //we really have no idea
	return match;
}
