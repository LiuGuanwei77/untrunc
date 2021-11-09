#include "codec.h"
#include "log.h"

extern "C" {
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#ifdef WIN32

#define H264_MAX_PICTURE_COUNT 36

#define MAX_MMCO_COUNT         66
#define MAX_DELAYED_PIC_COUNT  16
/**
 * The maximum number of slices supported by the decoder.
 * must be a power of 2
 */
#define MAX_SLICES 32
#define MAX_SPS_COUNT          32
#define MAX_PPS_COUNT         256

 /**
  * Sequence parameter set
  */
	typedef struct SPS {
		unsigned int sps_id;
		int profile_idc;
		int level_idc;
		int chroma_format_idc;
		int transform_bypass;              ///< qpprime_y_zero_transform_bypass_flag
		int log2_max_frame_num;            ///< log2_max_frame_num_minus4 + 4
		int poc_type;                      ///< pic_order_cnt_type
		int log2_max_poc_lsb;              ///< log2_max_pic_order_cnt_lsb_minus4
		int delta_pic_order_always_zero_flag;
		int offset_for_non_ref_pic;
		int offset_for_top_to_bottom_field;
		int poc_cycle_length;              ///< num_ref_frames_in_pic_order_cnt_cycle
		int ref_frame_count;               ///< num_ref_frames
		int gaps_in_frame_num_allowed_flag;
		int mb_width;                      ///< pic_width_in_mbs_minus1 + 1
		///< (pic_height_in_map_units_minus1 + 1) * (2 - frame_mbs_only_flag)
		int mb_height;
		int frame_mbs_only_flag;
		int mb_aff;                        ///< mb_adaptive_frame_field_flag
		int direct_8x8_inference_flag;
		int crop;                          ///< frame_cropping_flag

		/* those 4 are already in luma samples */
		unsigned int crop_left;            ///< frame_cropping_rect_left_offset
		unsigned int crop_right;           ///< frame_cropping_rect_right_offset
		unsigned int crop_top;             ///< frame_cropping_rect_top_offset
		unsigned int crop_bottom;          ///< frame_cropping_rect_bottom_offset
		int vui_parameters_present_flag;
		AVRational sar;
		int video_signal_type_present_flag;
		int full_range;
		int colour_description_present_flag;
		enum AVColorPrimaries color_primaries;
		enum AVColorTransferCharacteristic color_trc;
		enum AVColorSpace colorspace;
		int timing_info_present_flag;
		uint32_t num_units_in_tick;
		uint32_t time_scale;
		int fixed_frame_rate_flag;
		short offset_for_ref_frame[256]; // FIXME dyn aloc?
		int bitstream_restriction_flag;
		int num_reorder_frames;
		int scaling_matrix_present;
		uint8_t scaling_matrix4[6][16];
		uint8_t scaling_matrix8[6][64];
		int nal_hrd_parameters_present_flag;
		int vcl_hrd_parameters_present_flag;
		int pic_struct_present_flag;
		int time_offset_length;
		int cpb_cnt;                          ///< See H.264 E.1.2
		int initial_cpb_removal_delay_length; ///< initial_cpb_removal_delay_length_minus1 + 1
		int cpb_removal_delay_length;         ///< cpb_removal_delay_length_minus1 + 1
		int dpb_output_delay_length;          ///< dpb_output_delay_length_minus1 + 1
		int bit_depth_luma;                   ///< bit_depth_luma_minus8 + 8
		int bit_depth_chroma;                 ///< bit_depth_chroma_minus8 + 8
		int residual_color_transform_flag;    ///< residual_colour_transform_flag
		int constraint_set_flags;             ///< constraint_set[0-3]_flag
		uint8_t data[4096];
		size_t data_size;
	} SPS;

	typedef struct VideoDSPContext {
		/**
		 * Copy a rectangular area of samples to a temporary buffer and replicate
		 * the border samples.
		 *
		 * @param dst destination buffer
		 * @param dst_stride number of bytes between 2 vertically adjacent samples
		 *                   in destination buffer
		 * @param src source buffer
		 * @param dst_linesize number of bytes between 2 vertically adjacent
		 *                     samples in the destination buffer
		 * @param src_linesize number of bytes between 2 vertically adjacent
		 *                     samples in both the source buffer
		 * @param block_w width of block
		 * @param block_h height of block
		 * @param src_x x coordinate of the top left sample of the block in the
		 *                source buffer
		 * @param src_y y coordinate of the top left sample of the block in the
		 *                source buffer
		 * @param w width of the source buffer
		 * @param h height of the source buffer
		 */
		void (*emulated_edge_mc)(uint8_t* dst, const uint8_t* src,
			ptrdiff_t dst_linesize,
			ptrdiff_t src_linesize,
			int block_w, int block_h,
			int src_x, int src_y, int w, int h);

		/**
		 * Prefetch memory into cache (if supported by hardware).
		 *
		 * @param buf    pointer to buffer to prefetch memory from
		 * @param stride distance between two lines of buf (in bytes)
		 * @param h      number of lines to prefetch
		 */
		void (*prefetch)(uint8_t* buf, ptrdiff_t stride, int h);
	} VideoDSPContext;

	/**
	 * Context for storing H.264 DSP functions
	 */
	typedef struct H264DSPContext {
		/* weighted MC */
		void* weight_h264_pixels_tab[4];
		void* biweight_h264_pixels_tab[4];

		/* loop filter */
		void (*h264_v_loop_filter_luma)(uint8_t* pix /*align 16*/, int stride,
			int alpha, int beta, int8_t* tc0);
		void (*h264_h_loop_filter_luma)(uint8_t* pix /*align 4 */, int stride,
			int alpha, int beta, int8_t* tc0);
		void (*h264_h_loop_filter_luma_mbaff)(uint8_t* pix /*align 16*/, int stride,
			int alpha, int beta, int8_t* tc0);
		/* v/h_loop_filter_luma_intra: align 16 */
		void (*h264_v_loop_filter_luma_intra)(uint8_t* pix, int stride,
			int alpha, int beta);
		void (*h264_h_loop_filter_luma_intra)(uint8_t* pix, int stride,
			int alpha, int beta);
		void (*h264_h_loop_filter_luma_mbaff_intra)(uint8_t* pix /*align 16*/,
			int stride, int alpha, int beta);
		void (*h264_v_loop_filter_chroma)(uint8_t* pix /*align 8*/, int stride,
			int alpha, int beta, int8_t* tc0);
		void (*h264_h_loop_filter_chroma)(uint8_t* pix /*align 4*/, int stride,
			int alpha, int beta, int8_t* tc0);
		void (*h264_h_loop_filter_chroma_mbaff)(uint8_t* pix /*align 8*/,
			int stride, int alpha, int beta,
			int8_t* tc0);
		void (*h264_v_loop_filter_chroma_intra)(uint8_t* pix /*align 8*/,
			int stride, int alpha, int beta);
		void (*h264_h_loop_filter_chroma_intra)(uint8_t* pix /*align 8*/,
			int stride, int alpha, int beta);
		void (*h264_h_loop_filter_chroma_mbaff_intra)(uint8_t* pix /*align 8*/,
			int stride, int alpha, int beta);
		// h264_loop_filter_strength: simd only. the C version is inlined in h264_loopfilter.c
		void (*h264_loop_filter_strength)(int16_t bS[2][4][4], uint8_t nnz[40],
			int8_t ref[2][40], int16_t mv[2][40][2],
			int bidir, int edges, int step,
			int mask_mv0, int mask_mv1, int field);

		/* IDCT */
		void (*h264_idct_add)(uint8_t* dst /*align 4*/,
			int16_t* block /*align 16*/, int stride);
		void (*h264_idct8_add)(uint8_t* dst /*align 8*/,
			int16_t* block /*align 16*/, int stride);
		void (*h264_idct_dc_add)(uint8_t* dst /*align 4*/,
			int16_t* block /*align 16*/, int stride);
		void (*h264_idct8_dc_add)(uint8_t* dst /*align 8*/,
			int16_t* block /*align 16*/, int stride);

		void (*h264_idct_add16)(uint8_t* dst /*align 16*/, const int* blockoffset,
			int16_t* block /*align 16*/, int stride,
			const uint8_t nnzc[15 * 8]);
		void (*h264_idct8_add4)(uint8_t* dst /*align 16*/, const int* blockoffset,
			int16_t* block /*align 16*/, int stride,
			const uint8_t nnzc[15 * 8]);
		void (*h264_idct_add8)(uint8_t** dst /*align 16*/, const int* blockoffset,
			int16_t* block /*align 16*/, int stride,
			const uint8_t nnzc[15 * 8]);
		void (*h264_idct_add16intra)(uint8_t* dst /*align 16*/, const int* blockoffset,
			int16_t* block /*align 16*/,
			int stride, const uint8_t nnzc[15 * 8]);
		void (*h264_luma_dc_dequant_idct)(int16_t* output,
			int16_t* input /*align 16*/, int qmul);
		void (*h264_chroma_dc_dequant_idct)(int16_t* block, int qmul);

		/* bypass-transform */
		void (*h264_add_pixels8_clear)(uint8_t* dst, int16_t* block, int stride);
		void (*h264_add_pixels4_clear)(uint8_t* dst, int16_t* block, int stride);

		/**
		 * Search buf from the start for up to size bytes. Return the index
		 * of a zero byte, or >= size if not found. Ideally, use lookahead
		 * to filter out any zero bytes that are known to not be followed by
		 * one or more further zero bytes and a one byte. Better still, filter
		 * out any bytes that form the trailing_zero_8bits syntax element too.
		 */
		int (*startcode_find_candidate)(const uint8_t* buf, int size);
	} H264DSPContext;

	typedef void (*h264_chroma_mc_func)(uint8_t* dst /*align 8*/, uint8_t* src /*align 1*/, ptrdiff_t srcStride, int h, int x, int y);

	typedef struct H264ChromaContext {
		h264_chroma_mc_func put_h264_chroma_pixels_tab[4];
		h264_chroma_mc_func avg_h264_chroma_pixels_tab[4];
	} H264ChromaContext;

	typedef struct H264QpelContext {
		void* put_h264_qpel_pixels_tab[4][16];
		void* avg_h264_qpel_pixels_tab[4][16];
	} H264QpelContext;

	typedef struct ThreadFrame {
		AVFrame* f;
		AVCodecContext* owner[2];
		// progress->data is an array of 2 ints holding progress for top/bottom
		// fields
		AVBufferRef* progress;
	} ThreadFrame;

	typedef struct H264Picture {
		AVFrame* f;
		ThreadFrame tf;

		AVBufferRef* qscale_table_buf;
		int8_t* qscale_table;

		AVBufferRef* motion_val_buf[2];
		int16_t(*motion_val[2])[2];

		AVBufferRef* mb_type_buf;
		uint32_t* mb_type;

		AVBufferRef* hwaccel_priv_buf;
		void* hwaccel_picture_private; ///< hardware accelerator private data

		AVBufferRef* ref_index_buf[2];
		int8_t* ref_index[2];

		int field_poc[2];       ///< top/bottom POC
		int poc;                ///< frame POC
		int frame_num;          ///< frame_num (raw frame_num from slice header)
		int mmco_reset;         /**< MMCO_RESET set this 1. Reordering code must
									 not mix pictures before and after MMCO_RESET. */
		int pic_id;             /**< pic_num (short -> no wrap version of pic_num,
									 pic_num & max_pic_num; long -> long_pic_num) */
		int long_ref;           ///< 1->long term reference 0->short term reference
		int ref_poc[2][2][32];  ///< POCs of the frames/fields used as reference (FIXME need per slice)
		int ref_count[2][2];    ///< number of entries in ref_poc         (FIXME need per slice)
		int mbaff;              ///< 1 -> MBAFF frame 0-> not MBAFF
		int field_picture;      ///< whether or not picture was encoded in separate fields

		int reference;
		int recovered;          ///< picture at IDR or recovery point + recovery count
		int invalid_gap;
		int sei_recovery_frame_cnt;
	} H264Picture;

	struct H264SliceContext;
	struct H2645NAL;

	typedef struct H2645RBSP {
		uint8_t* rbsp_buffer;
		int rbsp_buffer_alloc_size;
		int rbsp_buffer_size;
	} H2645RBSP;

	/* an input packet split into unescaped NAL units */
	typedef struct H2645Packet {
		H2645NAL* nals;
		H2645RBSP rbsp;
		int nb_nals;
		int nals_allocated;
	} H2645Packet;

	/**
	 * Context for storing H.264 prediction functions
	 */
	typedef struct H264PredContext {
		void(*pred4x4[9 + 3 + 3])(uint8_t* src, const uint8_t* topright,
			ptrdiff_t stride);
		void(*pred8x8l[9 + 3])(uint8_t* src, int topleft, int topright,
			ptrdiff_t stride);
		void(*pred8x8[4 + 3 + 4])(uint8_t* src, ptrdiff_t stride);
		void(*pred16x16[4 + 3 + 2])(uint8_t* src, ptrdiff_t stride);

		void(*pred4x4_add[2])(uint8_t* pix /*align  4*/,
			int16_t* block /*align 16*/, ptrdiff_t stride);
		void(*pred8x8l_add[2])(uint8_t* pix /*align  8*/,
			int16_t* block /*align 16*/, ptrdiff_t stride);
		void(*pred8x8l_filter_add[2])(uint8_t* pix /*align  8*/,
			int16_t* block /*align 16*/, int topleft, int topright, ptrdiff_t stride);
		void(*pred8x8_add[3])(uint8_t* pix /*align  8*/,
			const int* block_offset,
			int16_t* block /*align 16*/, ptrdiff_t stride);
		void(*pred16x16_add[3])(uint8_t* pix /*align 16*/,
			const int* block_offset,
			int16_t* block /*align 16*/, ptrdiff_t stride);
	} H264PredContext;

	struct PPS;

	typedef struct H264ParamSets {
		AVBufferRef* sps_list[MAX_SPS_COUNT];
		AVBufferRef* pps_list[MAX_PPS_COUNT];

		AVBufferRef* pps_ref;
		AVBufferRef* sps_ref;
		/* currently active parameters sets */
		const PPS* pps;
		const SPS* sps;
	} H264ParamSets;

	typedef struct H264POCContext {
		int poc_lsb;
		int poc_msb;
		int delta_poc_bottom;
		int delta_poc[2];
		int frame_num;
		int prev_poc_msb;           ///< poc_msb of the last reference pic for POC type 0
		int prev_poc_lsb;           ///< poc_lsb of the last reference pic for POC type 0
		int frame_num_offset;       ///< for POC type 2
		int prev_frame_num_offset;  ///< for POC type 2
		int prev_frame_num;         ///< frame_num of the last pic for POC type 1/2
	} H264POCContext;

	typedef struct H264Ref {
		uint8_t* data[3];
		int linesize[3];

		int reference;
		int poc;
		int pic_id;

		H264Picture* parent;
	} H264Ref;

	/**
	 * Memory management control operation opcode.
	 */
	typedef enum MMCOOpcode {
		MMCO_END = 0,
		MMCO_SHORT2UNUSED,
		MMCO_LONG2UNUSED,
		MMCO_SHORT2LONG,
		MMCO_SET_MAX_LONG,
		MMCO_RESET,
		MMCO_LONG,
	} MMCOOpcode;

	/**
	 * Memory management control operation.
	 */
	typedef struct MMCO {
		MMCOOpcode opcode;
		int short_pic_num;  ///< pic_num without wrapping (pic_num & max_pic_num)
		int long_arg;       ///< index, pic_num, or num long refs depending on opcode
	} MMCO;

	/**
	 * H264Context
	 */
	typedef struct H264Context {
		const AVClass* _class;
		AVCodecContext* avctx;
		VideoDSPContext vdsp;
		H264DSPContext h264dsp;
		H264ChromaContext h264chroma;
		H264QpelContext h264qpel;

		H264Picture DPB[H264_MAX_PICTURE_COUNT];
		H264Picture* cur_pic_ptr;
		H264Picture cur_pic;
		H264Picture last_pic_for_ec;

		H264SliceContext* slice_ctx;
		int            nb_slice_ctx;
		int            nb_slice_ctx_queued;

		H2645Packet pkt;

		int pixel_shift;    ///< 0 for 8-bit H.264, 1 for high-bit-depth H.264

		/* coded dimensions -- 16 * mb w/h */
		int width, height;
		int chroma_x_shift, chroma_y_shift;

		int droppable;
		int coded_picture_number;

		int context_initialized;
		int flags;
		int workaround_bugs;
		int x264_build;
		/* Set when slice threading is used and at least one slice uses deblocking
		 * mode 1 (i.e. across slice boundaries). Then we disable the loop filter
		 * during normal MB decoding and execute it serially at the end.
		 */
		int postpone_filter;

		/*
		 * Set to 1 when the current picture is IDR, 0 otherwise.
		 */
		int picture_idr;

		int crop_left;
		int crop_right;
		int crop_top;
		int crop_bottom;

		int8_t(*intra4x4_pred_mode);
		H264PredContext hpc;

		uint8_t(*non_zero_count)[48];

#define LIST_NOT_USED -1 // FIXME rename?
#define PART_NOT_AVAILABLE -2

		/**
		 * block_offset[ 0..23] for frame macroblocks
		 * block_offset[24..47] for field macroblocks
		 */
		int block_offset[2 * (16 * 3)];

		uint32_t* mb2b_xy;  // FIXME are these 4 a good idea?
		uint32_t* mb2br_xy;
		int b_stride;       // FIXME use s->b4_stride

		uint16_t* slice_table;      ///< slice_table_base + 2*mb_stride + 1

		// interlacing specific flags
		int mb_aff_frame;
		int picture_structure;
		int first_field;

		uint8_t* list_counts;               ///< Array of list_count per MB specifying the slice type

		/* 0x100 -> non null luma_dc, 0x80/0x40 -> non null chroma_dc (cb/cr), 0x?0 -> chroma_cbp(0, 1, 2), 0x0? luma_cbp */
		uint16_t* cbp_table;

		/* chroma_pred_mode for i4x4 or i16x16, else 0 */
		uint8_t* chroma_pred_mode_table;
		uint8_t(*mvd_table[2])[2];
		uint8_t* direct_table;

		uint8_t scan_padding[16];
		uint8_t zigzag_scan[16];
		uint8_t zigzag_scan8x8[64];
		uint8_t zigzag_scan8x8_cavlc[64];
		uint8_t field_scan[16];
		uint8_t field_scan8x8[64];
		uint8_t field_scan8x8_cavlc[64];
		uint8_t zigzag_scan_q0[16];
		uint8_t zigzag_scan8x8_q0[64];
		uint8_t zigzag_scan8x8_cavlc_q0[64];
		uint8_t field_scan_q0[16];
		uint8_t field_scan8x8_q0[64];
		uint8_t field_scan8x8_cavlc_q0[64];

		int mb_y;
		int mb_height, mb_width;
		int mb_stride;
		int mb_num;

		// =============================================================
		// Things below are not used in the MB or more inner code

		int nal_ref_idc;
		int nal_unit_type;

		int has_slice;          ///< slice NAL is found in the packet, set by decode_nal_units, its state does not need to be preserved outside h264_decode_frame()

		/**
		 * Used to parse AVC variant of H.264
		 */
		int is_avc;           ///< this flag is != 0 if codec is avc1
		int nal_length_size;  ///< Number of bytes used for nal length (1, 2 or 4)

		int bit_depth_luma;         ///< luma bit depth from sps to detect changes
		int chroma_format_idc;      ///< chroma format from sps to detect changes

		H264ParamSets ps;

		uint16_t* slice_table_base;

		H264POCContext poc;

		H264Ref default_ref[2];
		H264Picture* short_ref[32];
		H264Picture* long_ref[32];
		H264Picture* delayed_pic[MAX_DELAYED_PIC_COUNT + 2]; // FIXME size?
		int last_pocs[MAX_DELAYED_PIC_COUNT];
		H264Picture* next_output_pic;
		int next_outputed_poc;

		/**
		 * memory management control operations buffer.
		 */
		MMCO mmco[MAX_MMCO_COUNT];
		int  nb_mmco;
		int mmco_reset;
		int explicit_ref_marking;

		int long_ref_count;     ///< number of actual long term references
		int short_ref_count;    ///< number of actual short term references

		/**
		 * @name Members for slice based multithreading
		 * @{
		 */
		 /**
		  * current slice number, used to initialize slice_num of each thread/context
		  */
		int current_slice;

		/** @} */

		/**
		 * Complement sei_pic_struct
		 * SEI_PIC_STRUCT_TOP_BOTTOM and SEI_PIC_STRUCT_BOTTOM_TOP indicate interlaced frames.
		 * However, soft telecined frames may have these values.
		 * This is used in an attempt to flag soft telecine progressive.
		 */
		int prev_interlaced_frame;

		/**
		 * Are the SEI recovery points looking valid.
		 */
		int valid_recovery_point;

		/**
		 * recovery_frame is the frame_num at which the next frame should
		 * be fully constructed.
		 *
		 * Set to -1 when not expecting a recovery point.
		 */
		int recovery_frame;

		/**
		 * We have seen an IDR, so all the following frames in coded order are correctly
		 * decodable.
		 */
#define FRAME_RECOVERED_IDR  (1 << 0)
		 /**
		  * Sufficient number of frames have been decoded since a SEI recovery point,
		  * so all the following frames in presentation order are correct.
		  */
#define FRAME_RECOVERED_SEI  (1 << 1)

		int frame_recovered;    ///< Initial frame has been completely recovered

		int has_recovery_point;

		int missing_fields;

		/* for frame threading, this is set to 1
		 * after finish_setup() has been called, so we cannot modify
		 * some context properties (which are supposed to stay constant between
		 * slices) anymore */
		int setup_finished;

		int cur_chroma_format_idc;
		int cur_bit_depth_luma;
		int16_t slice_row[MAX_SLICES]; ///< to detect when MAX_SLICES is too low

		/* original AVCodecContext dimensions, used to handle container
		 * cropping */
		int width_from_caller;
		int height_from_caller;

		int enable_er;

		//H264SEIContext sei;

		AVBufferPool* qscale_table_pool;
		AVBufferPool* mb_type_pool;
		AVBufferPool* motion_val_pool;
		AVBufferPool* ref_index_pool;
		int ref2frm[MAX_SLICES][2][64];     ///< reference to frame number lists, used in the loop filter, the first 2 are for -2,-1
	} H264Context;
#else
#include <libavcodec/h264dec.h>
#endif
}
#include <assert.h>

using namespace std;




// AVC1
class H264sps {
public:
	int  log2_max_frame_num;
	bool frame_mbs_only_flag;
	int  poc_type;
	int  log2_max_poc_lsb;

	H264sps()
		: log2_max_frame_num(0),
		  frame_mbs_only_flag(false),
		  poc_type(0),
		  log2_max_poc_lsb(0)
	{ }

	H264sps(const SPS &avsps)
		: log2_max_frame_num(avsps.log2_max_frame_num),
		  frame_mbs_only_flag(bool(avsps.frame_mbs_only_flag)),
		  poc_type(avsps.poc_type),
		  log2_max_poc_lsb(avsps.log2_max_poc_lsb)
	{ }

	H264sps &operator=(const SPS &avsps) { return *this = H264sps(avsps); }

	void parseSPS(const uint8_t *data, int size);
};


void H264sps::parseSPS(const uint8_t *data, int size) {

	assert(data != NULL);
	if (data[0] != 1) {
		Log::debug << "Uncharted territory...\n";
	}

	if (size < 7) {
		throw string("Could not parse SPS!");
	}

	// Decode SPS from avcC.
	const uint8_t *p   = data;
	int            cnt = p[5] & 0x1f;   // Number of SPS.
	p += 6;
	if(cnt != 1) {
		Log::debug << "Not supporting more than 1 SPS unit for the moment; might fail horribly.\n";
	}
	for(int i = 0; i < cnt; i++) {
		//nalsize = AV_RB16(p) + 2;
		int nalsize = readBE<uint16_t>(p);
		if (p - data + nalsize > size) {
			throw string("Could not parse SPS!");
		}
#if 0
		// From: libavcodec/h264_parse.c
		ret = decode_extradata_ps_mp4(p, nalsize, ps, err_recognition, logctx);
		if (ret < 0) {
			av_log(logctx, AV_LOG_ERROR,
				   "Decoding SPS %d from avcC failed\n", i);
			return ret;
		}
		p += nalsize;
#endif
		break;
	}

	// Skip PPS.
}


class NalInfo {

public:
	static const int MaxAVC1Length = 8 * (1 << 20);

	int length;

	int ref_idc;
	int nal_type;
	int first_mb;           // Unused.
	int slice_type;         // Should match the NAL type (1, 5).
	int pps_id;             // Pic parameter set id: which parameter set to use.
	int frame_num;
	int field_pic_flag;
	int bottom_pic_flag;
	int idr_pic_flag;       // Actually 1 for nal_type 5, 0 for nal_type 0.
	int idr_pic_id;         // Read only for nal_type 5.
	int poc_type;           // If zero, check the poc lsb.
	int poc_lsb;            // Pic order count - least significant bit.

	NalInfo()
		: length(0),
		  ref_idc(0),
		  nal_type(0),
		  first_mb(0),
		  slice_type(0),
		  pps_id(0),
		  frame_num(0),
		  field_pic_flag(0),
		  bottom_pic_flag(0),
		  idr_pic_flag(0),
		  idr_pic_id(0),
		  poc_type(0),
		  poc_lsb(0)
	{ }

	bool getNalInfo(const H264sps &sps, uint32_t maxlength, const uint8_t *buffer);
	void clear();
	void print(int indentation = 0);

private:
	static int golomb  (       uint8_t *&buffer, int &offset);
	static int readBits(int n, uint8_t *&buffer, int &offset);
};


void NalInfo::clear() {
	length          = 0;
	ref_idc         = 0;
	nal_type        = 0;
	first_mb        = 0;
	slice_type      = 0;
	pps_id          = 0;
	frame_num       = 0;
	field_pic_flag  = 0;
	bottom_pic_flag = 0;
	idr_pic_flag    = 0;
	idr_pic_id      = 0;
	poc_type        = 0;
	poc_lsb         = 0;
}

void NalInfo::print(int indentation) {
	const string indent((indentation >= 0)? indentation : 0, ' ');

	Log::debug << indent << "Length         : " << length          << ((length < 8+4 || length > MaxAVC1Length+4) ? " (incorrect)":"") << '\n';
	Log::debug << indent << "Ref idc        : " << ref_idc         << ((unsigned(ref_idc) > 3) ? " (incorrect)":"") << '\n';
	Log::debug << indent << "Nal type       : " << nal_type        << ((unsigned(nal_type) > 0x1f) ? " (incorrect)":"") << '\n';
	Log::debug << indent << "First mb       : " << first_mb        << '\n';
	Log::debug << indent << "Slice type     : " << slice_type      << ((unsigned(slice_type) > 9) ? " (incorrect)":"") << '\n';
	Log::debug << indent << "Pic parm set id: " << pps_id          << '\n';
	Log::debug << indent << "Frame number   : " << frame_num       << '\n';
	Log::debug << indent << "Field  pic flag: " << field_pic_flag  << '\n';
	if(field_pic_flag)
		Log::debug << indent << "Bottom pic flag: " << bottom_pic_flag << '\n';
	Log::debug << indent << "Idr pic id     : " << idr_pic_id      << '\n';
	if(poc_type)
		Log::debug << indent << "Poc type       : " << poc_type    << '\n';
	else
		Log::debug << indent << "Poc lsb        : " << poc_lsb     << '\n';
	Log::flush();
}


int NalInfo::golomb(uint8_t *&buffer, int &offset) {
	assert(buffer != NULL && offset >= 0 && offset < 8);
	// Count the leading zeroes.
	int count = 0;
	while((*buffer & (0x01 << (7 - offset))) == 0) {
		count++;
		offset++;
		if(offset == 8) {
			buffer++;
			offset = 0;
		}
		if(count > 20) {
			cerr << "Failed reading golomb: too large!\n";
			return -1;
		}
	}
	// Skip the single 1 delimiter.
	offset++;
	if(offset == 8) {
		buffer++;
		offset = 0;
	}
	uint32_t res = 1;
	// Read count bits.
	while(count-- > 0) {
		res <<= 1;
		res |= (*buffer & (0x01 << (7 - offset))) >> (7 - offset);
		//res |= (*buffer >> (7 - offset)) & 0x01;
		offset++;
		if(offset == 8) {
			buffer++;
			offset = 0;
		}
	}
	return res - 1;
}

int NalInfo::readBits(int n, uint8_t *&buffer, int &offset) {
	assert(buffer != NULL && offset >= 0);
	int res = 0;
	// Can't read in a single reading.
	while(n + offset > 8) {
		int d = 8 - offset;
		res <<= d;
		res |= *buffer & ((1 << d) - 1);
		offset = 0;
		buffer++;
		n -= d;
	}
	// Read the remaining bits.
	int d = (8 - offset - n);
	res <<= n;
	res |= (*buffer >> d) & ((1 << n) - 1);
	return res;
}

// Return false means this probably is not a NAL.
bool NalInfo::getNalInfo(const H264sps &sps, uint32_t maxlength, const uint8_t *buffer) {
	// Re-initialize.
	clear();

	if(buffer[0] != 0) {
		Log::debug << "First byte expected 0.\n";
		return false;
	}

	// This is supposed to be the length of the NAL unit.
	uint32_t len = readBE<uint32_t>(buffer);
	if(len > MaxAVC1Length) {
		Log::debug << "Max AVC1 length exceeded.\n";
		return false;
	}
	if(len > maxlength) {
		Log::debug << "Buffer size exceeded (" << (len ) << " > " << maxlength << ").\n";
		return false;
	}
	length = len + 4;
	Log::debug << "Length         : " << length << '\n';

	buffer += 4;
	if(*buffer & (1 << 7)) {
		Log::debug << "Forbidden first bit 1.\n";
		return false;
	}
	ref_idc = *buffer >> 5;
	//Log::debug << "Ref idc        : " << ref_idc << '\n';

	nal_type = *buffer & 0x1f;
	Log::debug << "Nal type       : " << nal_type << '\n';
	if(nal_type != 1 && nal_type != 5)
		return true;

	// Check if size is reasonable.
	if(len < 6) {
		Log::debug << "Length too short! (" << len << " < 7).\n";
		return false;
	}

	// Skip NAL header.
	buffer++;

	// Remove the emulation prevention 0x03 byte.
	// Could be done in place to speed up things.
	vector<uint8_t> data;
	data.reserve(len);
	for(unsigned int i = 0; i < len; i++) {
		data.push_back(buffer[i]);
		if(i+2 < len && buffer[i] == 0 && buffer[i+1] == 0 && buffer[i+2] == 3) {
			data.push_back(buffer[i+1]);
			assert(buffer[i+2] == 0x03);
			i += 2; // Skipping 3 byte!
		}
	}

	uint8_t *start  = data.data();
	int      offset = 0;
	first_mb   = golomb(start, offset);
	if(first_mb < 0) return false;
	// TODO: Is there a max number, so we could validate?
	//Log::debug << "First mb       : " << first_mb << '\n';

	slice_type = golomb(start, offset);
	if(slice_type < 0) return false;

	if(slice_type > 9) {
		Log::debug << "Invalid slice type (" << slice_type << "), probably this is not an avc1 sample.\n";
		return false;
	}
	Log::debug << "Slice type     : " << slice_type << '\n';

	pps_id     = golomb(start, offset);
	if(pps_id < 0) return false;

	//Log::debug << "Pic parm set id: " << pps_id << '\n';
	// pps id: should be taked from master context (h264_slice.c:1257).

	// Assume separate colour plane flag is 0,
	//  otherwise we would have to read colour_plane_id which is 2 bits.

	// Assuming same sps for all frames.
	//SPS *sps = reinterpret_cast<SPS *>(h->ps.sps_list[0]->data);  // may_alias.
	frame_num = readBits(sps.log2_max_frame_num, start, offset);
	Log::debug << "Frame number   : " << frame_num << '\n';

	// Read 2 flags.
	field_pic_flag  = 0;
	bottom_pic_flag = 0;
	if(sps.frame_mbs_only_flag) {
		field_pic_flag = readBits(1, start, offset);
		//Log::debug << "Field  pic flag: " << field_pic_flag << '\n';
		if(field_pic_flag) {
			bottom_pic_flag = readBits(1, start, offset);
			//Log::debug << "Bottom pic flag: " << bottom_pic_flag << '\n';
		}
	}

	idr_pic_flag = (nal_type == 5) ? 1 : 0;
	if (idr_pic_flag) {
		idr_pic_id = golomb(start, offset);
		if(idr_pic_id < 0) return false;

		//Log::debug << "Idr pic id     : " << idr_pic_id << '\n';
	}

	// If the pic order count type == 0.
	poc_type = sps.poc_type;
	if(sps.poc_type == 0) {
		poc_lsb = readBits(sps.log2_max_poc_lsb, start, offset);
		//Log::debug << "Poc lsb        : " << poc_lsb << '\n';
	}

	// Ignoring the delta_poc for the moment.
	return true;
}





Match Codec::avc1Search(const unsigned char *start, int maxlength) {
	Match match;
	for(int offset = 0; offset < maxlength - 8; offset++) {
		if(start[offset] != 0)
			continue;
		uint32_t len = readBE<uint32_t>(start + offset);
		//too might use smallestSample and largestSample to constrain the size
		if(len < 8 || len > NalInfo::MaxAVC1Length)
			continue;
		//todo common values for 4 and 5 bytes should be usedc
		//if(start[offset+4] != 0x41 || start[offset+5] != 0x9a ) continue;
		if(start[offset + 4] & (1 << 7))
			continue;
		int nal_type = start[offset + 4] & 0x1f;
		if(nal_type > 21)
			continue;
		//this looks like it might be a packet. might want to do a more thorough check.
		match.offset = offset;
		match.chances = 10;
		return match;
	}
	return match;
}


Match Codec::avc1Match(const unsigned char *start, int maxlength) {
	Match match;

	// This works only for a very specific kind of video.
	//#define SPECIAL_VIDEO
#ifdef SPECIAL_VIDEO
	int32_t s2 = readBE<int32_t>(start + 4);
	if(s != 0x00000002 || (s2 != 0x09300000 && s2 != 0x09100000))
		return false;
	return true;
#endif

	//First 4 bytes is the length the the packet and it's expected not to be bigger than 16M
	if(start[0] != 0) {
//		Log::debug << "avc1: Match with 0 header.\n";
		return match;
	}

	// NAL unit types
	//  see: libavcodec/h264.h
	//  see: ITU-T T-REC-H.264-201704, Table 7-1
	// enum {
	//     NAL_SLICE             =  1,  // Non keyframe.
	//     NAL_DPA               =  2,
	//     NAL_DPB               =  3,
	//     NAL_DPC               =  4,
	//     NAL_IDR_SLICE         =  5,  // Keyframe.
	//     NAL_SEI               =  6,
	//     NAL_SPS               =  7,
	//     NAL_PPS               =  8,
	//     NAL_AUD               =  9,
	//     NAL_END_SEQUENCE      = 10,
	//     NAL_END_STREAM        = 11,
	//     NAL_FILLER_DATA       = 12,
	//     NAL_SPS_EXT           = 13,
	//     NAL_PREFIX            = 14,
	//     NAL_SUB_SPS           = 15,
	//     NAL_DPS               = 16,
	//     NAL_AUXILIARY_SLICE   = 19,best
	//     NAL_EXTEN_SLICE       = 20,
	//     NAL_DEPTH_EXTEN_SLICE = 21,
	//
	//     NAL_FF_IGNORE         = 0xff0f001,
	// };
	//
	// First 4 bytes are the length, then the NAL starts.
	// ref_idc != 0 per unit_type = 5
	// ref_idc == 0 per unit_type = 6, 9, 10, 11, 12

	// See 7.4.1.2.4 Detection of the first VCL NAL unit of a primary coded picture
	//  for rules on how to group NALs into a picture.

	int64_t begin32 = readBE<int32_t>(start);
	if(begin32 == 0)
		return match;

	// TODO: Use the first byte of the NAL: forbidden bit and type!
	int nal_type = (start[4] & 0x1f);

	if(nal_type == 0)
		return match;
	// The other values are really uncommon on cameras...
	if(nal_type > 12) {
		Log::debug << "avc1: No match because of NAL type: " << nal_type << '\n';
		return match;
	}
	if(nal_type > 29) {
		//if(nal_type != 1 && !(nal_type >= 5 && nal_type <= 12)) {
		Log::debug << "avc1: No match because of NAL type: " << nal_type << '\n';
		return match;
	}

	if(!context) {
		Log::error << "avc1: No context" << nal_type << '\n';
        return match;
	}

	// XXX: Horrible Hack: Referencing unstable, internal data structures. XXX
	H264Context *h    = static_cast<H264Context*>(context->priv_data); //context->codec->
	const SPS   *hsps = NULL;
	if(h) {
		// Use currently active SPS.
		if(!hsps && h->ps.sps_list[0]) {
			// Use first SPS.
			hsps = reinterpret_cast<const SPS*>(h->ps.sps_list[0]->data); // may_alias.
		}
	}
	if(!hsps) {
		Log::debug << "Could not retrieve SPS.\n";
		return match;
	}
	H264sps sps(*hsps);

#if 0
	int consumed = -1;
	{
		AvLog useAvLog();
		AVFrame *frame = av_frame_alloc();
		if(!frame)
			throw string("Could not create AVFrame");
		AVPacket avp;
		av_init_packet(&avp);
		avp.data = start;
		avp.size = maxlength;
		int got_frame = 0;
		consumed = avcodec_decode_video2(context, frame, &got_frame, &avp);
		if(consumed == 0) {
			// Flush decoder to receive buffered packets.
			clog << "Flush " << name << " decoder.\n";
			got_frame = 0;
			av_packet_unref(&avp);
			av_frame_unref(frame);
			int consumed2 = avcodec_decode_video2(context, frame, &got_frame, &avp);
			if(consumed2 >= 0)
				consumed = consumed2;
		}
		av_packet_unref(&avp);
		av_frame_free(&frame);
	}

	return consumed;
#endif


	uint32_t  length = 0;
	const unsigned char *pos    = start;

	NalInfo previous;
	bool  seen_slice = false;

	bool first_pack = true;
	while(true) {
		NalInfo info;
		bool ok = info.getNalInfo(sps, maxlength, pos);
		if(!ok) {
			//THIS should never happens, but it happens
			if(first_pack) {

				//throw string("What's ahppinghin egherklhj HELP!");
				match.chances = 0.0f;
				if(info.length == 0) {
					NalInfo info1;

					info1.getNalInfo(sps, maxlength, pos);
					return match;
				}
				return match;
			}
			goto final;
		}
		match.chances = 1024;

		first_pack = false;

		Log::debug << "Nal type: " << info.nal_type << " Nal ref idc " << info.ref_idc << endl;

		switch(info.nal_type) {
		case 12: break; //filler data
		case 1:
		case 5:
			if(!seen_slice) {
				previous   = info;
				seen_slice = true;
			} else {
				// Check for changes.
				//cout << "Frame number: " << info.frame_num << endl;
				if(previous.frame_num != info.frame_num) {
					Log::debug << "Different frame number.\n";
					goto final;

				}
				if(previous.pps_id != info.pps_id) {
					Log::debug << "Different pic parameter set id.\n";
					goto final;
				}
				// All these conditions are listed in the docs, but
				// it looks like it creates invalid packets if respected.
				// Puzzling.
				//#define STRICT_NAL_INFO_CHECKING  1
//#define STRICT_FIELD_PIC_CHECKING
#ifdef STRICT_FIELD_PIC_CHECKING
				if(previous.field_pic_flag != info.field_pic_flag) {
					Log::debug << "Different field  pic flag.\n";
					goto final;
				}
				if(previous.field_pic_flag && info.field_pic_flag
						&& previous.bottom_pic_flag != info.bottom_pic_flag)
				{
					Log::debug << "Different bottom pic flag.\n";
					goto final;
				}
#endif

#define STRICT_REF_IDC_CHECKING
#ifdef STRICT_REF_IDC_CHECKING
				if(previous.ref_idc != info.ref_idc) {
					Log::debug << "Different ref idc.\n";
					goto final;
				}
#endif

//#define STRICT_POC_TYPE_CHECKING
#ifdef STRICT_POC_TYPE_CHECKING
				if(previous.poc_type == 0 && info.poc_type == 0
						&& previous.poc_lsb != info.poc_lsb)
				{
					Log::debug << "Different pic order count lsb (poc lsb).\n";
					goto final;
				}
#endif
				if(previous.idr_pic_flag != info.idr_pic_flag) {
					Log::debug << "Different NAL type (5, 1).\n";
					//TODO check this was an error and not on purpouse.
					goto final;
				}
//#define STRICT_PIC_IDR_CHECKING
#ifdef STRICT_PIC_IDR_CHECKING
				if(previous.idr_pic_flag == 1 && info.idr_pic_flag == 1
						&& previous.idr_pic_id != info.idr_pic_id)
				{
					Log::debug << "Different idr pic id for keyframe.\n";
					goto final;
				}
#endif
			}
			break;
		default:
			if(seen_slice) {
				Log::debug << "New access unit since seen picture.\n";
				goto final;
			}
			break;
		}
		if(info.nal_type == 5) {
			match.keyframe = true;
			Log::debug << "KeyFrame!" << endl;
		}
		pos       += info.length;
		length    += info.length;
		maxlength -= info.length;
		//Log::debug << "Partial length : " << length << '\n';
	}
	final:
	match.length = length;
	//probability depends on length, the longer the higher.
	match.chances = 1 + length/10;
	if(maxlength < 8)
		return match;

	float stat_chanches = 0;
	if(stats.beginnings32.count(begin32))
		stat_chanches = stats.beginnings32[begin32];
	else {
		//TODO this actually depends by the number of different beginnings.
		//changes = (n -1); //se sono due le chanches sono davvero piccole.
		stat_chanches = stats.beginnings32.size() -1;
	}

	int64_t begin64 = readBE<int64_t>(start);
	if(stats.beginnings64.count(begin64))
		stat_chanches = stats.beginnings64[begin64];

	match.chances = std::max(stat_chanches, match.chances);

	Log::flush();
	return match;
}
