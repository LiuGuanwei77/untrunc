#include "codec.h"
#include "codec_handler.h"
#include "log.h"
#include <string.h>

using namespace std;

/* mbex is a codec for some editing data I cannot find documentation for.
 * for the few examples I got they are very small packets (8, 10, 100 bytes
 * large packets included crec and cist atoms.
 * Examples:
 *
 * Length 8: 00 00 00 08   00 00 00 00   01 3E 15 29
 * Length 10: 00 00 00 0A   00 00 00 01   00 06 00 00 23 95
 * Length 100: 00 00 00 64    00 00 00 01   00 00 00 5C  63 72 65 63  00 00 00 54  63 69 74 73  00 00 00 0C  00 00 00 02  00 00 00 01
 */

namespace {
// Refactor step 28: bodies inlined from the now-removed Codec::apchMatch
// and Codec::apchSearch. No Codec members are used; the codec ref is
// unused.
struct ApchHandler : CodecHandler {
	Match match(Codec &, const unsigned char *start, int maxlength) override {
		Match match;
		uint32_t length = readBE<uint32_t>(start);

		if(!strncmp((char *)start + 4, "icpf", 4)) {
			match.length = length;
			match.chances  = 1e30;
			return match;
		}
		return match;
	}
	Match search(Codec &, const unsigned char *start, int maxlength) override {
		Match match;
		const unsigned char *end = start + maxlength - 8;
		const  unsigned char *current = start;
		while(current < end) {
			if(!strncmp((char *)current + 4, "icpf", 4)) {
				match.offset = current - start;
				match.chances  = 1e30;
				return match;
			}
			current++;
		}
		return match;
	}
};
} // namespace

std::unique_ptr<CodecHandler> makeApchHandler() {
	return std::unique_ptr<CodecHandler>(new ApchHandler());
}
