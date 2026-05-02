#ifndef CODEC_HANDLER_H
#define CODEC_HANDLER_H

class Codec;
struct Match;

// Abstract base for per-codec match/search strategies.
// During the in-progress refactor, concrete subclasses delegate back to the
// existing Codec::xxxMatch / Codec::xxxSearch member functions so behavior is
// byte-identical to the pre-refactor dispatcher. Once every codec has a
// handler, the bodies will be inlined here and the Codec::xxxMatch members
// removed.
class CodecHandler {
public:
	virtual ~CodecHandler() = default;
	virtual Match match(Codec &codec, const unsigned char *start, int maxlength) = 0;
	virtual Match search(Codec &codec, const unsigned char *start, int maxlength);
};

#endif // CODEC_HANDLER_H
