#include "codec_handler.h"
#include "codec.h"

// Default search returns an empty Match — preserves the pre-refactor behavior
// of Codec::search() falling through for codecs that don't override it.
Match CodecHandler::search(Codec &, const unsigned char *, int) {
	return Match();
}
