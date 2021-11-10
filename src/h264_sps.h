#ifndef H264_SPS_H
#define H264_SPS_H


struct H264sps {
  int log2_max_frame_num;
  int frame_mbs_only_flag;
  int poc_type;
  int log2_max_poc_lsb;
};

int initH264sps(struct H264sps* sps, const void* context);
  
#endif