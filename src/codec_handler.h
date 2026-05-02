#ifndef CODEC_HANDLER_H
#define CODEC_HANDLER_H

#include <memory>

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

// Per-codec handler factories. During the migration, each handler delegates
// back to the existing Codec::xxxMatch / Codec::xxxSearch member function so
// dispatch is byte-identical to the pre-refactor behavior.
std::unique_ptr<CodecHandler> makeApchHandler();
std::unique_ptr<CodecHandler> makeCammHandler();
std::unique_ptr<CodecHandler> makeFdscHandler();
std::unique_ptr<CodecHandler> makeGpmdHandler();
std::unique_ptr<CodecHandler> makeMbexHandler();
std::unique_ptr<CodecHandler> makeMijdHandler();
std::unique_ptr<CodecHandler> makeRtpHandler();
std::unique_ptr<CodecHandler> makeTextHandler();
std::unique_ptr<CodecHandler> makeTmcdHandler();
std::unique_ptr<CodecHandler> makeTx3gHandler();

#endif // CODEC_HANDLER_H
