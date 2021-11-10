#include "h264_sps.h"

#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>

#include "libavcodec/h264dec.h"

void SPS2H264sps(struct H264sps* sps, const SPS* avsps) {
  sps->log2_max_frame_num = avsps->log2_max_frame_num;
  sps->frame_mbs_only_flag = avsps->frame_mbs_only_flag;
  sps->poc_type = avsps->poc_type;
  sps->log2_max_poc_lsb = avsps->log2_max_poc_lsb;
}

int initH264sps(struct H264sps* sps, const void* context) {
  H264Context* h = (H264Context*)context;
  const SPS* hsps = NULL;
  if (h) {
    // Use currently active SPS.
    if (!hsps && h->ps.sps_list[0]) {
      // Use first SPS.
      hsps = (SPS*)h->ps.sps_list[0]->data;
    }
  }
  if (!hsps) {
    printf("Could not retrieve SPS.\n");
    return FALSE;
  }
  SPS2H264sps(sps, hsps);
  return TRUE;
}
