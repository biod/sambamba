module cram.htslib;

// TODO: structs with bitfields are not represented correctly

import core.stdc.stdio;
import core.stdc.config;
import core.stdc.stdarg;
import std.c.linux.pthread;

extern (C):

alias __kstring_t kstring_t;
alias _Anonymous_0 htsFile;
alias __hts_idx_t hts_idx_t;
alias _Anonymous_1 hts_pair64_t;
alias int function (BGZF*, void*, void*, int*, int*, int*) hts_readrec_func;
alias _Anonymous_2 hts_itr_t;
alias int function (void*, const(char)*) hts_name2id_f;
alias _Anonymous_2* function (const(__hts_idx_t)*, int, int, int, int function (BGZF*, void*, void*, int*, int*, int*)) hts_itr_query_func;
alias const(char)* function (void*, int) hts_id2name_f;
alias _Anonymous_3 bam_hdr_t;
alias _Anonymous_4 bam1_core_t;
alias _Anonymous_5 bam1_t;
alias _Anonymous_0 samFile;
alias _Anonymous_6 bam_pileup1_t;
alias int function (void*, _Anonymous_5*) bam_plp_auto_f;
alias __bam_plp_t* bam_plp_t;
alias __bam_mplp_t* bam_mplp_t;
alias _Anonymous_5 bam_seq_t;
alias _Anonymous_7 string_t;
alias _Anonymous_8 string_alloc_t;
alias _Anonymous_9 pool_t;
alias _Anonymous_10 pool_alloc_t;
alias uint khint32_t;
alias c_ulong khint64_t;
alias uint khint_t;
alias uint khiter_t;
alias const(char)* kh_cstr_t;
alias _Anonymous_11 ks_tokaux_t;
alias SAM_hdr_tag_s SAM_hdr_tag;
alias SAM_hdr_item_s SAM_hdr_type;
alias _Anonymous_12 SAM_SQ;
alias _Anonymous_13 SAM_RG;
alias _Anonymous_14 SAM_PG;
alias _Anonymous_15 kh_sam_hdr_t;
alias _Anonymous_16 kh_m_s2i_t;
alias _Anonymous_17 SAM_hdr;
alias t_res t_pool_result;
alias _Anonymous_18 t_pool_worker_t;
alias _Anonymous_19 kh_m_i2i_t;
alias _Anonymous_20 kh_s_i2i_t;
alias ubyte uc;
alias _Anonymous_21 pmap_t;
alias _Anonymous_22 kh_map_t;
alias _Anonymous_23 cram_stats;
alias _Anonymous_24 cram_file_def;
alias _Anonymous_25 cram_metrics;
alias _Anonymous_26 cram_block;
alias _Anonymous_27 cram_block_compression_hdr;
alias _Anonymous_28 cram_block_slice_hdr;
alias _Anonymous_29 cram_container;
alias _Anonymous_30 cram_record;
alias _Anonymous_31 cram_feature;
alias _Anonymous_32 kh_refs_t;
alias _Anonymous_33 refs_t;
alias _Anonymous_34 cram_range;
alias _Anonymous_35 cram_huffman_code;
alias _Anonymous_36 cram_huffman_decoder;
alias _Anonymous_37 cram_huffman_encoder;
alias _Anonymous_38 cram_beta_decoder;
alias _Anonymous_39 cram_gamma_decoder;
alias _Anonymous_40 cram_subexp_decoder;
alias _Anonymous_41 cram_external_decoder;
alias _Anonymous_42 cram_byte_array_len_decoder;
alias _Anonymous_43 cram_byte_array_stop_decoder;
alias _Anonymous_44 cram_byte_array_len_encoder;

extern __gshared int hts_verbose;
extern __gshared const ubyte[256] seq_nt16_table;
extern __gshared const ubyte* seq_nt16_str;
extern __gshared const double __ac_HASH_UPPER;

enum htsFormatCategory
{
	unknown_category = 0,
	sequence_data = 1,
	variant_data = 2,
	index_file = 3,
	region_list = 4,
	category_maximum = 32767
}

enum htsExactFormat
{
	unknown_format = 0,
	binary_format = 1,
	text_format = 2,
	sam = 3,
	bam = 4,
	bai = 5,
	cram = 6,
	crai = 7,
	vcf = 8,
	bcfv1 = 9,
	bcf = 10,
	csi = 11,
	gzi = 12,
	tbi = 13,
	bed = 14,
	format_maximum = 32767
}

enum htsCompression
{
	no_compression = 0,
	gzip = 1,
	bgzf = 2,
	custom = 3,
	compression_maximum = 32767
}

enum sam_fields
{
	SAM_QNAME = 1,
	SAM_FLAG = 2,
	SAM_RNAME = 4,
	SAM_POS = 8,
	SAM_MAPQ = 16,
	SAM_CIGAR = 32,
	SAM_RNEXT = 64,
	SAM_PNEXT = 128,
	SAM_TLEN = 256,
	SAM_SEQ = 512,
	SAM_QUAL = 1024,
	SAM_AUX = 2048,
	SAM_RGAUX = 4096
}

enum cram_option
{
	CRAM_OPT_DECODE_MD = 0,
	CRAM_OPT_PREFIX = 1,
	CRAM_OPT_VERBOSITY = 2,
	CRAM_OPT_SEQS_PER_SLICE = 3,
	CRAM_OPT_SLICES_PER_CONTAINER = 4,
	CRAM_OPT_RANGE = 5,
	CRAM_OPT_VERSION = 6,
	CRAM_OPT_EMBED_REF = 7,
	CRAM_OPT_IGNORE_MD5 = 8,
	CRAM_OPT_REFERENCE = 9,
	CRAM_OPT_MULTI_SEQ_PER_SLICE = 10,
	CRAM_OPT_NO_REF = 11,
	CRAM_OPT_USE_BZIP2 = 12,
	CRAM_OPT_SHARED_REF = 13,
	CRAM_OPT_NTHREADS = 14,
	CRAM_OPT_THREAD_POOL = 15,
	CRAM_OPT_USE_LZMA = 16,
	CRAM_OPT_USE_RANS = 17,
	CRAM_OPT_REQUIRED_FIELDS = 18
}

enum cigar_op
{
	BAM_CMATCH_ = 0,
	BAM_CINS_ = 1,
	BAM_CDEL_ = 2,
	BAM_CREF_SKIP_ = 3,
	BAM_CSOFT_CLIP_ = 4,
	BAM_CHARD_CLIP_ = 5,
	BAM_CPAD_ = 6,
	BAM_CBASE_MATCH = 7,
	BAM_CBASE_MISMATCH = 8
}

enum cram_encoding
{
	E_NULL = 0,
	E_EXTERNAL = 1,
	E_GOLOMB = 2,
	E_HUFFMAN = 3,
	E_BYTE_ARRAY_LEN = 4,
	E_BYTE_ARRAY_STOP = 5,
	E_BETA = 6,
	E_SUBEXP = 7,
	E_GOLOMB_RICE = 8,
	E_GAMMA = 9
}

enum cram_external_type
{
	E_INT = 1,
	E_LONG = 2,
	E_BYTE = 3,
	E_BYTE_ARRAY = 4,
	E_BYTE_ARRAY_BLOCK = 5
}

enum cram_DS_ID
{
	DS_CORE = 0,
	DS_aux = 1,
	DS_aux_OQ = 2,
	DS_aux_BQ = 3,
	DS_aux_BD = 4,
	DS_aux_BI = 5,
	DS_aux_FZ = 6,
	DS_aux_oq = 7,
	DS_aux_os = 8,
	DS_aux_oz = 9,
	DS_ref = 10,
	DS_RN = 11,
	DS_QS = 12,
	DS_IN = 13,
	DS_SC = 14,
	DS_BF = 15,
	DS_CF = 16,
	DS_AP = 17,
	DS_RG = 18,
	DS_MQ = 19,
	DS_NS = 20,
	DS_MF = 21,
	DS_TS = 22,
	DS_NP = 23,
	DS_NF = 24,
	DS_RL = 25,
	DS_FN = 26,
	DS_FC = 27,
	DS_FP = 28,
	DS_DL = 29,
	DS_BA = 30,
	DS_BS = 31,
	DS_TL = 32,
	DS_RI = 33,
	DS_RS = 34,
	DS_PD = 35,
	DS_HC = 36,
	DS_BB = 37,
	DS_QQ = 38,
	DS_TN = 39,
	DS_RN_len = 40,
	DS_SC_len = 41,
	DS_BB_len = 42,
	DS_QQ_len = 43,
	DS_TC = 44,
	DS_TM = 45,
	DS_TV = 46,
	DS_END = 47
}

enum cram_block_method
{
	ERROR = -1,
	RAW = 0,
	GZIP = 1,
	BZIP2 = 2,
	LZMA = 3,
	RANS = 4,
	RANS0 = 4,
	RANS1 = 10,
	GZIP_RLE = 11
}

enum cram_content_type
{
	CT_ERROR = -1,
	FILE_HEADER = 0,
	COMPRESSION_HEADER = 1,
	MAPPED_SLICE = 2,
	UNMAPPED_SLICE = 3,
	EXTERNAL = 4,
	CORE = 5
}

enum cram_fields
{
	CRAM_BF = 1,
	CRAM_AP = 2,
	CRAM_FP = 4,
	CRAM_RL = 8,
	CRAM_DL = 16,
	CRAM_NF = 32,
	CRAM_BA = 64,
	CRAM_QS = 128,
	CRAM_FC = 256,
	CRAM_FN = 512,
	CRAM_BS = 1024,
	CRAM_IN = 2048,
	CRAM_RG = 4096,
	CRAM_MQ = 8192,
	CRAM_TL = 16384,
	CRAM_RN = 32768,
	CRAM_NS = 65536,
	CRAM_NP = 131072,
	CRAM_TS = 262144,
	CRAM_MF = 524288,
	CRAM_CF = 1048576,
	CRAM_RI = 2097152,
	CRAM_RS = 4194304,
	CRAM_PD = 8388608,
	CRAM_HC = 16777216,
	CRAM_SC = 33554432,
	CRAM_BB = 67108864,
	CRAM_BB_len = 134217728,
	CRAM_QQ = 268435456,
	CRAM_QQ_len = 536870912,
	CRAM_aux = 1073741824,
	CRAM_ALL = 2147483647
}

struct __kstring_t
{
	size_t l;
	size_t m;
	char* s;
}

struct htsFormat
{
	enum htsFormatCategory
	{
		unknown_category = 0,
		sequence_data = 1,
		variant_data = 2,
		index_file = 3,
		region_list = 4,
		category_maximum = 32767
	}
	htsFormatCategory category;
	enum htsExactFormat
	{
		unknown_format = 0,
		binary_format = 1,
		text_format = 2,
		sam = 3,
		bam = 4,
		bai = 5,
		cram = 6,
		crai = 7,
		vcf = 8,
		bcfv1 = 9,
		bcf = 10,
		csi = 11,
		gzi = 12,
		tbi = 13,
		bed = 14,
		format_maximum = 32767
	}
	htsExactFormat format;
	enum htsCompression
	{
		no_compression = 0,
		gzip = 1,
		bgzf = 2,
		custom = 3,
		compression_maximum = 32767
	}
	htsCompression compression;
}

struct _Anonymous_0
{
	uint is_bin;
	uint is_write;
	uint is_be;
	uint is_cram;
	uint dummy;
	long lineno;
	kstring_t line;
	char* fn;
	char* fn_aux;
	union
	{
		BGZF* bgzf;
		cram_fd* cram;
		hFILE* hfile;
		void* voidp;
	}
	htsFormat format;
}

struct _Anonymous_1
{
	ulong u;
	ulong v;
}

struct _Anonymous_2
{
	uint read_rest;
	uint finished;
	uint dummy;
	int tid;
	int beg;
	int end;
	int n_off;
	int i;
	ulong curr_off;
	hts_pair64_t* off;
	int function (BGZF*, void*, void*, int*, int*, int*) readrec;
	struct
	{
		int n;
		int m;
		int* a;
	}
}

struct _Anonymous_3
{
	int n_targets;
	int ignore_sam_err;
	uint l_text;
	uint* target_len;
	byte* cigar_tab;
	char** target_name;
	char* text;
	void* sdict;
}

struct _Anonymous_4
{
	int tid;
	int pos;
	uint bin_mq_nl;
	uint flag_nc;
	int l_qseq;
	int mtid;
	int mpos;
	int isize;
}

struct _Anonymous_5
{
	bam1_core_t core;
	int l_data;
	int m_data;
	ubyte* data;
	ulong id;
}

struct _Anonymous_6
{
	bam1_t* b;
	int qpos;
	int indel;
	int level;
	uint is_del;
	uint is_head;
	uint is_tail;
	uint is_refskip;
	uint aux;
}

struct _Anonymous_7
{
	char* str;
	size_t used;
}

struct _Anonymous_8
{
	size_t max_length;
	size_t nstrings;
	string_t* strings;
}

struct _Anonymous_9
{
	void* pool;
	size_t used;
}

struct _Anonymous_10
{
	size_t dsize;
	size_t npools;
	pool_t* pools;
	void* free;
}

struct _Anonymous_11
{
	ulong[4] tab;
	int sep;
	int finished;
	const(char)* p;
}

struct SAM_hdr_tag_s
{
	SAM_hdr_tag_s* next;
	char* str;
	int len;
}

struct SAM_hdr_item_s
{
	SAM_hdr_item_s* next;
	SAM_hdr_item_s* prev;
	SAM_hdr_tag* tag;
	int order;
}

struct _Anonymous_12
{
	char* name;
	uint len;
	SAM_hdr_type* ty;
	SAM_hdr_tag* tag;
}

struct _Anonymous_13
{
	char* name;
	SAM_hdr_type* ty;
	SAM_hdr_tag* tag;
	int name_len;
	int id;
}

struct _Anonymous_14
{
	char* name;
	SAM_hdr_type* ty;
	SAM_hdr_tag* tag;
	int name_len;
	int id;
	int prev_id;
}

struct _Anonymous_15
{
	khint_t n_buckets;
	khint_t size;
	khint_t n_occupied;
	khint_t upper_bound;
	khint32_t* flags;
	khint32_t* keys;
	SAM_hdr_type** vals;
}

struct _Anonymous_16
{
	khint_t n_buckets;
	khint_t size;
	khint_t n_occupied;
	khint_t upper_bound;
	khint32_t* flags;
	kh_cstr_t* keys;
	int* vals;
}

struct _Anonymous_17
{
	kstring_t text;
	kh_sam_hdr_t* h;
	string_alloc_t* str_pool;
	pool_alloc_t* type_pool;
	pool_alloc_t* tag_pool;
	int nref;
	SAM_SQ* ref_;
	kh_m_s2i_t* ref_hash;
	int nrg;
	SAM_RG* rg;
	kh_m_s2i_t* rg_hash;
	int npg;
	int npg_end;
	int npg_end_alloc;
	SAM_PG* pg;
	kh_m_s2i_t* pg_hash;
	int* pg_end;
	char[1024] ID_buf;
	int ID_cnt;
	int ref_count;
}

struct t_pool_job
{
	void* function (void*) func;
	void* arg;
	t_pool_job* next;
	t_pool* p;
	t_results_queue* q;
	int serial;
}

struct t_res
{
	t_res* next;
	int serial;
	void* data;
}

struct _Anonymous_18
{
	t_pool* p;
	int idx;
	pthread_t tid;
	pthread_cond_t pending_c;
	long wait_time;
}

struct t_pool
{
	int qsize;
	int njobs;
	int nwaiting;
	int shutdown;
	t_pool_job* head;
	t_pool_job* tail;
	int tsize;
	t_pool_worker_t* t;
	pthread_mutex_t pool_m;
	pthread_cond_t empty_c;
	pthread_cond_t pending_c;
	pthread_cond_t full_c;
	int* t_stack;
	int t_stack_top;
	long total_time;
	long wait_time;
}

struct t_results_queue
{
	t_pool_result* result_head;
	t_pool_result* result_tail;
	int next_serial;
	int curr_serial;
	int queue_len;
	int pending;
	pthread_mutex_t result_m;
	pthread_cond_t result_avail_c;
}

struct _Anonymous_19
{
	khint_t n_buckets;
	khint_t size;
	khint_t n_occupied;
	khint_t upper_bound;
	khint32_t* flags;
	khint32_t* keys;
	int* vals;
}

struct _Anonymous_20
{
	khint_t n_buckets;
	khint_t size;
	khint_t n_occupied;
	khint_t upper_bound;
	khint32_t* flags;
	khint32_t* keys;
	char* vals;
}

struct _Anonymous_22
{
	khint_t n_buckets;
	khint_t size;
	khint_t n_occupied;
	khint_t upper_bound;
	khint32_t* flags;
	kh_cstr_t* keys;
	pmap_t* vals;
}

struct _Anonymous_23
{
	int[1024] freqs;
	kh_m_i2i_t* h;
	int nsamp;
	int nvals;
}

struct _Anonymous_24
{
	char[4] magic;
	ubyte major_version;
	ubyte minor_version;
	char[20] file_id;
}

struct _Anonymous_25
{
	int trial;
	int next_trial;
	int sz_gz_rle;
	int sz_gz_def;
	int sz_rans0;
	int sz_rans1;
	int sz_bzip2;
	int sz_lzma;
	int method;
	int strat;
	int gz_rle_cnt;
	int gz_def_cnt;
	int rans0_cnt;
	int rans1_cnt;
	int bzip2_cnt;
	int lzma_cnt;
	int revised_method;
	double gz_rle_extra;
	double gz_def_extra;
	double rans0_extra;
	double rans1_extra;
	double bzip2_extra;
	double lzma_extra;
}

struct _Anonymous_26
{
	enum cram_block_method
	{
		ERROR = -1,
		RAW = 0,
		GZIP = 1,
		BZIP2 = 2,
		LZMA = 3,
		RANS = 4,
		RANS0 = 4,
		RANS1 = 10,
		GZIP_RLE = 11
	}
	cram_block_method method;
	cram_block_method orig_method;
	enum cram_content_type
	{
		CT_ERROR = -1,
		FILE_HEADER = 0,
		COMPRESSION_HEADER = 1,
		MAPPED_SLICE = 2,
		UNMAPPED_SLICE = 3,
		EXTERNAL = 4,
		CORE = 5
	}
	cram_content_type content_type;
	int content_id;
	int comp_size;
	int uncomp_size;
	uint crc32;
	int idx;
	ubyte* data;
	size_t alloc;
	size_t byte_;
	int bit;
}

struct _Anonymous_27
{
	int ref_seq_id;
	int ref_seq_start;
	int ref_seq_span;
	int num_records;
	int num_landmarks;
	int* landmark;
	int mapped_qs_included;
	int unmapped_qs_included;
	int unmapped_placed;
	int qs_included;
	int read_names_included;
	int AP_delta;
	char[4][5] substitution_matrix;
	cram_block* TD_blk;
	int nTL;
	ubyte** TL;
	kh_m_s2i_t* TD_hash;
	string_alloc_t* TD_keys;
	kh_map_t* preservation_map;
	cram_map*[32] rec_encoding_map;
	cram_map*[32] tag_encoding_map;
	cram_codec*[47] codecs;
	char* uncomp;
	size_t uncomp_size;
	size_t uncomp_alloc;
	uint data_series;
}

struct cram_map
{
	int key;
	enum cram_encoding
	{
		E_NULL = 0,
		E_EXTERNAL = 1,
		E_GOLOMB = 2,
		E_HUFFMAN = 3,
		E_BYTE_ARRAY_LEN = 4,
		E_BYTE_ARRAY_STOP = 5,
		E_BETA = 6,
		E_SUBEXP = 7,
		E_GOLOMB_RICE = 8,
		E_GAMMA = 9
	}
	cram_encoding encoding;
	int offset;
	int size;
	cram_codec* codec;
	cram_map* next;
}

struct _Anonymous_28
{
	enum cram_content_type
	{
		CT_ERROR = -1,
		FILE_HEADER = 0,
		COMPRESSION_HEADER = 1,
		MAPPED_SLICE = 2,
		UNMAPPED_SLICE = 3,
		EXTERNAL = 4,
		CORE = 5
	}
	cram_content_type content_type;
	int ref_seq_id;
	int ref_seq_start;
	int ref_seq_span;
	int num_records;
	int record_counter;
	int num_blocks;
	int num_content_ids;
	int* block_content_ids;
	int ref_base_id;
	ubyte[16] md5;
}

struct _Anonymous_29
{
	int length;
	int ref_seq_id;
	int ref_seq_start;
	int ref_seq_span;
	int record_counter;
	long num_bases;
	int num_records;
	int num_blocks;
	int num_landmarks;
	int* landmark;
	size_t offset;
	cram_block_compression_hdr* comp_hdr;
	cram_block* comp_hdr_block;
	int max_slice;
	int curr_slice;
	int max_rec;
	int curr_rec;
	int max_c_rec;
	int curr_c_rec;
	int slice_rec;
	int curr_ref;
	int last_pos;
	cram_slice** slices;
	cram_slice* slice;
	int pos_sorted;
	int max_apos;
	int last_slice;
	int multi_seq;
	int unsorted;
	int ref_start;
	int first_base;
	int last_base;
	int ref_id;
	int ref_end;
	char* ref_;
	bam_seq_t** bams;
	cram_stats*[47] stats;
	kh_s_i2i_t* tags_used;
	int* refs_used;
	uint crc32;
}

struct _Anonymous_30
{
	cram_slice* s;
	int ref_id;
	int flags;
	int cram_flags;
	int len;
	int apos;
	int rg;
	int name;
	int name_len;
	int mate_line;
	int mate_ref_id;
	int mate_pos;
	int tlen;
	int ntags;
	int aux;
	int aux_size;
	int TN_idx;
	int TL;
	int seq;
	int qual;
	int cigar;
	int ncigar;
	int aend;
	int mqual;
	int feature;
	int nfeature;
	int mate_flags;
}

struct _Anonymous_31
{
}

struct cram_slice
{
	cram_block_slice_hdr* hdr;
	cram_block* hdr_block;
	cram_block** block;
	cram_block** block_by_id;
	int last_apos;
	int max_apos;
	ulong id;
	cram_record* crecs;
	uint* cigar;
	uint cigar_alloc;
	uint ncigar;
	cram_feature* features;
	int nfeatures;
	int afeatures;
	uint* TN;
	int nTN;
	int aTN;
	cram_block* name_blk;
	cram_block* seqs_blk;
	cram_block* qual_blk;
	cram_block* base_blk;
	cram_block* soft_blk;
	cram_block* aux_blk;
	cram_block* aux_OQ_blk;
	cram_block* aux_BQ_blk;
	cram_block* aux_BD_blk;
	cram_block* aux_BI_blk;
	cram_block* aux_FZ_blk;
	cram_block* aux_oq_blk;
	cram_block* aux_os_blk;
	cram_block* aux_oz_blk;
	string_alloc_t* pair_keys;
	kh_m_s2i_t*[2] pair;
	char* ref_;
	int ref_start;
	int ref_end;
	int ref_id;
}

struct ref_entry
{
	char* name;
	char* fn;
	long length;
	long offset;
	int bases_per_line;
	int line_length;
	long count;
	char* seq;
}

struct _Anonymous_32
{
	khint_t n_buckets;
	khint_t size;
	khint_t n_occupied;
	khint_t upper_bound;
	khint32_t* flags;
	kh_cstr_t* keys;
	ref_entry** vals;
}

struct _Anonymous_33
{
	string_alloc_t* pool;
	kh_refs_t* h_meta;
	ref_entry** ref_id;
	int nref;
	char* fn;
	BGZF* fp;
	int count;
	pthread_mutex_t lock;
	ref_entry* last;
	int last_id;
}

struct cram_index
{
	int nslice;
	int nalloc;
	cram_index* e;
	int refid;
	int start;
	int end;
	int nseq;
	int slice;
	int len;
	long offset;
}

struct _Anonymous_34
{
	int refid;
	int start;
	int end;
}

struct spare_bams
{
	bam_seq_t** bams;
	spare_bams* next;
}

struct cram_fd
{
	hFILE* fp;
	int mode;
	int version_;
	cram_file_def* file_def;
	SAM_hdr* header;
	char* prefix;
	int record_counter;
	int slice_num;
	int err;
	cram_container* ctr;
	int first_base;
	int last_base;
	refs_t* refs;
	char* ref_;
	char* ref_free;
	int ref_id;
	int ref_start;
	int ref_end;
	char* ref_fn;
	int level;
	cram_metrics*[47] m;
	int decode_md;
	int verbose;
	int seqs_per_slice;
	int slices_per_container;
	int embed_ref;
	int no_ref;
	int ignore_md5;
	int use_bz2;
	int use_rans;
	int use_lzma;
	int shared_ref;
	uint required_fields;
	cram_range range;
	uint[4096] bam_flag_swap;
	uint[4096] cram_flag_swap;
	ubyte[256] L1;
	ubyte[256] L2;
	char[32][32] cram_sub_matrix;
	int index_sz;
	cram_index* index;
	off_t first_container;
	int eof;
	int last_slice;
	int multi_seq;
	int unsorted;
	int empty_container;
	int own_pool;
	t_pool* pool;
	t_results_queue* rqueue;
	pthread_mutex_t metrics_lock;
	pthread_mutex_t ref_lock;
	spare_bams* bl;
	pthread_mutex_t bam_list_lock;
	void* job_pending;
	int ooc;
}

struct _Anonymous_35
{
	int symbol;
	int p;
	int code;
	int len;
}

struct _Anonymous_36
{
	int ncodes;
	cram_huffman_code* codes;
}

struct _Anonymous_37
{
	cram_huffman_code* codes;
	int nvals;
	int[129] val2code;
}

struct _Anonymous_38
{
	int offset;
	int nbits;
}

struct _Anonymous_39
{
	int offset;
}

struct _Anonymous_40
{
	int offset;
	int k;
}

struct _Anonymous_41
{
	int content_id;
	enum cram_external_type
	{
		E_INT = 1,
		E_LONG = 2,
		E_BYTE = 3,
		E_BYTE_ARRAY = 4,
		E_BYTE_ARRAY_BLOCK = 5
	}
	cram_external_type type;
}

struct _Anonymous_42
{
	cram_codec* len_codec;
	cram_codec* value_codec;
}

struct _Anonymous_43
{
	ubyte stop;
	int content_id;
}

struct _Anonymous_44
{
	enum cram_encoding
	{
		E_NULL = 0,
		E_EXTERNAL = 1,
		E_GOLOMB = 2,
		E_HUFFMAN = 3,
		E_BYTE_ARRAY_LEN = 4,
		E_BYTE_ARRAY_STOP = 5,
		E_BETA = 6,
		E_SUBEXP = 7,
		E_GOLOMB_RICE = 8,
		E_GAMMA = 9
	}
	cram_encoding len_encoding;
	cram_encoding val_encoding;
	void* len_dat;
	void* val_dat;
	cram_codec* len_codec;
	cram_codec* val_codec;
}

struct cram_codec
{
	enum cram_encoding
	{
		E_NULL = 0,
		E_EXTERNAL = 1,
		E_GOLOMB = 2,
		E_HUFFMAN = 3,
		E_BYTE_ARRAY_LEN = 4,
		E_BYTE_ARRAY_STOP = 5,
		E_BETA = 6,
		E_SUBEXP = 7,
		E_GOLOMB_RICE = 8,
		E_GAMMA = 9
	}
	cram_encoding codec;
	cram_block* out_;
	void function (cram_codec*) free;
	int function (cram_slice*, cram_codec*, cram_block*, char*, int*) decode;
	int function (cram_slice*, cram_codec*, char*, int) encode;
	int function (cram_codec*, cram_block*, char*, int) store;
}

struct hFILE;


struct __bam_plp_t;


struct BGZF;


struct __bam_mplp_t;


struct __hts_idx_t;


union _Anonymous_21
{
	int i;
	char* p;
}

const(char)* hts_version ();
int hts_detect_format (hFILE* fp, htsFormat* fmt);
const(char)* hts_format_description (const(htsFormat)* format);
htsFile* hts_open (const(char)* fn, const(char)* mode);
htsFile* hts_hopen (hFILE* fp, const(char)* fn, const(char)* mode);
int hts_close (htsFile* fp);
const(htsFormat)* hts_get_format (htsFile* fp);
int hts_set_opt (htsFile* fp, cram_option opt, ...);
int hts_getline (htsFile* fp, int delimiter, kstring_t* str);
char** hts_readlines (const(char)* fn, int* _n);
char** hts_readlist (const(char)* fn, int is_file, int* _n);
int hts_set_threads (htsFile* fp, int n);
int hts_set_fai_filename (htsFile* fp, const(char)* fn_aux);
hts_idx_t* hts_idx_init (int n, int fmt, ulong offset0, int min_shift, int n_lvls);
void hts_idx_destroy (hts_idx_t* idx);
int hts_idx_push (hts_idx_t* idx, int tid, int beg, int end, ulong offset, int is_mapped);
void hts_idx_finish (hts_idx_t* idx, ulong final_offset);
void hts_idx_save (const(hts_idx_t)* idx, const(char)* fn, int fmt);
hts_idx_t* hts_idx_load (const(char)* fn, int fmt);
ubyte* hts_idx_get_meta (hts_idx_t* idx, int* l_meta);
void hts_idx_set_meta (hts_idx_t* idx, int l_meta, ubyte* meta, int is_copy);
int hts_idx_get_stat (const(hts_idx_t)* idx, int tid, ulong* mapped, ulong* unmapped);
ulong hts_idx_get_n_no_coor (const(hts_idx_t)* idx);
const(char)* hts_parse_reg (const(char)* s, int* beg, int* end);
hts_itr_t* hts_itr_query (const(hts_idx_t)* idx, int tid, int beg, int end, int function (BGZF*, void*, void*, int*, int*, int*) readrec);
void hts_itr_destroy (hts_itr_t* iter);
hts_itr_t* hts_itr_querys (const(hts_idx_t)* idx, const(char)* reg, hts_name2id_f getid, void* hdr, hts_itr_t* function (const(hts_idx_t)*, int, int, int, int function (BGZF*, void*, void*, int*, int*, int*)) itr_query, int function (BGZF*, void*, void*, int*, int*, int*) readrec);
int hts_itr_next (BGZF* fp, hts_itr_t* iter, void* r, void* data);
const(char*)* hts_idx_seqnames (const(hts_idx_t)* idx, int* n, hts_id2name_f getid, void* hdr);
int hts_reg2bin (long beg, long end, int min_shift, int n_lvls);
int hts_bin_bot (int bin, int n_lvls);
int ed_is_big ();
ushort ed_swap_2 (ushort v);
void* ed_swap_2p (void* x);
uint ed_swap_4 (uint v);
void* ed_swap_4p (void* x);
ulong ed_swap_8 (ulong v);
void* ed_swap_8p (void* x);
bam_hdr_t* bam_hdr_init ();
bam_hdr_t* bam_hdr_read (BGZF* fp);
int bam_hdr_write (BGZF* fp, const(bam_hdr_t)* h);
void bam_hdr_destroy (bam_hdr_t* h);
int bam_name2id (bam_hdr_t* h, const(char)* ref_);
bam_hdr_t* bam_hdr_dup (const(bam_hdr_t)* h0);
bam1_t* bam_init1 ();
void bam_destroy1 (bam1_t* b);
int bam_read1 (BGZF* fp, bam1_t* b);
int bam_write1 (BGZF* fp, const(bam1_t)* b);
bam1_t* bam_copy1 (bam1_t* bdst, const(bam1_t)* bsrc);
bam1_t* bam_dup1 (const(bam1_t)* bsrc);
int bam_cigar2qlen (int n_cigar, const(uint)* cigar);
int bam_cigar2rlen (int n_cigar, const(uint)* cigar);
int bam_endpos (const(bam1_t)* b);
int bam_str2flag (const(char)* str);
char* bam_flag2str (int flag);
int bam_index_build (const(char)* fn, int min_shift);
hts_idx_t* sam_index_load (htsFile* fp, const(char)* fn);
hts_itr_t* sam_itr_queryi (const(hts_idx_t)* idx, int tid, int beg, int end);
hts_itr_t* sam_itr_querys (const(hts_idx_t)* idx, bam_hdr_t* hdr, const(char)* region);
int sam_open_mode (char* mode, const(char)* fn, const(char)* format);
bam_hdr_t* sam_hdr_parse (int l_text, const(char)* text);
bam_hdr_t* sam_hdr_read (samFile* fp);
int sam_hdr_write (samFile* fp, const(bam_hdr_t)* h);
int sam_parse1 (kstring_t* s, bam_hdr_t* h, bam1_t* b);
int sam_format1 (const(bam_hdr_t)* h, const(bam1_t)* b, kstring_t* str);
int sam_read1 (samFile* fp, bam_hdr_t* h, bam1_t* b);
int sam_write1 (samFile* fp, const(bam_hdr_t)* h, const(bam1_t)* b);
ubyte* bam_aux_get (const(bam1_t)* b, const char[2] tag);
int bam_aux2i (const(ubyte)* s);
double bam_aux2f (const(ubyte)* s);
char bam_aux2A (const(ubyte)* s);
char* bam_aux2Z (const(ubyte)* s);
void bam_aux_append (bam1_t* b, const char[2] tag, char type, int len, ubyte* data);
int bam_aux_del (bam1_t* b, ubyte* s);
bam_plp_t bam_plp_init (bam_plp_auto_f func, void* data);
void bam_plp_destroy (bam_plp_t iter);
int bam_plp_push (bam_plp_t iter, const(bam1_t)* b);
const(bam_pileup1_t)* bam_plp_next (bam_plp_t iter, int* _tid, int* _pos, int* _n_plp);
const(bam_pileup1_t)* bam_plp_auto (bam_plp_t iter, int* _tid, int* _pos, int* _n_plp);
void bam_plp_set_maxcnt (bam_plp_t iter, int maxcnt);
void bam_plp_reset (bam_plp_t iter);
bam_mplp_t bam_mplp_init (int n, bam_plp_auto_f func, void** data);
void bam_mplp_init_overlaps (bam_mplp_t iter);
void bam_mplp_destroy (bam_mplp_t iter);
void bam_mplp_set_maxcnt (bam_mplp_t iter, int maxcnt);
int bam_mplp_auto (bam_mplp_t iter, int* _tid, int* _pos, int* n_plp, const(bam_pileup1_t*)* plp);
string_alloc_t* string_pool_create (size_t max_length);
void string_pool_destroy (string_alloc_t* a_str);
char* string_alloc (string_alloc_t* a_str, size_t length);
char* string_dup (string_alloc_t* a_str, char* instr);
char* string_ndup (string_alloc_t* a_str, char* instr, size_t len);
pool_alloc_t* pool_create (size_t dsize);
void pool_destroy (pool_alloc_t* p);
void* pool_alloc (pool_alloc_t* p);
void pool_free (pool_alloc_t* p, void* ptr);
khint_t __ac_X31_hash_string (const(char)* s);
khint_t __ac_Wang_hash (khint_t key);
int kvsprintf (kstring_t* s, const(char)* fmt, va_list ap);
int ksprintf (kstring_t* s, const(char)* fmt, ...);
int ksplit_core (char* s, int delimiter, int* _max, int** _offsets);
char* kstrstr (const(char)* str, const(char)* pat, int** _prep);
char* kstrnstr (const(char)* str, const(char)* pat, int n, int** _prep);
void* kmemmem (const(void)* _str, int n, const(void)* _pat, int m, int** _prep);
char* kstrtok (const(char)* str, const(char)* sep, ks_tokaux_t* aux);
int ks_resize (kstring_t* s, size_t size);
char* ks_str (kstring_t* s);
size_t ks_len (kstring_t* s);
char* ks_release (kstring_t* s);
int kputsn (const(char)* p, int l, kstring_t* s);
int kputs (const(char)* p, kstring_t* s);
int kputc (int c, kstring_t* s);
int kputc_ (int c, kstring_t* s);
int kputsn_ (const(void)* p, int l, kstring_t* s);
int kputw (int c, kstring_t* s);
int kputuw (uint c, kstring_t* s);
int kputl (c_long c, kstring_t* s);
int* ksplit (kstring_t* s, int delimiter, int* n);
kh_sam_hdr_t* kh_init_sam_hdr ();
void kh_destroy_sam_hdr (kh_sam_hdr_t* h);
void kh_clear_sam_hdr (kh_sam_hdr_t* h);
khint_t kh_get_sam_hdr (const(kh_sam_hdr_t)* h, khint32_t key);
int kh_resize_sam_hdr (kh_sam_hdr_t* h, khint_t new_n_buckets);
khint_t kh_put_sam_hdr (kh_sam_hdr_t* h, khint32_t key, int* ret);
void kh_del_sam_hdr (kh_sam_hdr_t* h, khint_t x);
kh_m_s2i_t* kh_init_m_s2i ();
void kh_destroy_m_s2i (kh_m_s2i_t* h);
void kh_clear_m_s2i (kh_m_s2i_t* h);
khint_t kh_get_m_s2i (const(kh_m_s2i_t)* h, kh_cstr_t key);
int kh_resize_m_s2i (kh_m_s2i_t* h, khint_t new_n_buckets);
khint_t kh_put_m_s2i (kh_m_s2i_t* h, kh_cstr_t key, int* ret);
void kh_del_m_s2i (kh_m_s2i_t* h, khint_t x);
SAM_hdr* sam_hdr_new ();
SAM_hdr* sam_hdr_parse_ (const(char)* hdr, int len);
SAM_hdr* sam_hdr_dup (SAM_hdr* hdr);
void sam_hdr_incr_ref (SAM_hdr* hdr);
void sam_hdr_decr_ref (SAM_hdr* hdr);
void sam_hdr_free (SAM_hdr* hdr);
int sam_hdr_length (SAM_hdr* hdr);
char* sam_hdr_str (SAM_hdr* hdr);
int sam_hdr_add_lines (SAM_hdr* sh, const(char)* lines, int len);
int sam_hdr_add (SAM_hdr* sh, const(char)* type, ...);
int sam_hdr_vadd (SAM_hdr* sh, const(char)* type, va_list ap, ...);
SAM_hdr_type* sam_hdr_find (SAM_hdr* hdr, char* type, char* ID_key, char* ID_value);
char* sam_hdr_find_line (SAM_hdr* hdr, char* type, char* ID_key, char* ID_value);
SAM_hdr_tag* sam_hdr_find_key (SAM_hdr* sh, SAM_hdr_type* type, char* key, SAM_hdr_tag** prev);
int sam_hdr_update (SAM_hdr* hdr, SAM_hdr_type* type, ...);
int sam_hdr_rebuild (SAM_hdr* hdr);
int sam_hdr_name2ref (SAM_hdr* hdr, const(char)* ref_);
SAM_RG* sam_hdr_find_rg (SAM_hdr* hdr, const(char)* rg);
int sam_hdr_link_pg (SAM_hdr* hdr);
int sam_hdr_add_PG (SAM_hdr* sh, const(char)* name, ...);
bam_hdr_t* cram_header_to_bam (SAM_hdr* h);
SAM_hdr* bam_header_to_cram (bam_hdr_t* h);
int bam_construct_seq (bam_seq_t** bp, size_t extra_len, const(char)* qname, size_t qname_len, int flag, int rname, int pos, int end, int mapq, uint ncigar, const(uint)* cigar, int mrnm, int mpos, int isize, int len, const(char)* seq, const(char)* qual);
t_pool* t_pool_init (int qsize, int tsize);
int t_pool_dispatch (t_pool* p, t_results_queue* q, void* function (void*) func, void* arg);
int t_pool_dispatch2 (t_pool* p, t_results_queue* q, void* function (void*) func, void* arg, int nonblock);
int t_pool_flush (t_pool* p);
void t_pool_destroy (t_pool* p, int kill);
t_pool_result* t_pool_next_result (t_results_queue* q);
t_pool_result* t_pool_next_result_wait (t_results_queue* q);
void t_pool_delete_result (t_pool_result* r, int free_data);
t_results_queue* t_results_queue_init ();
void t_results_queue_destroy (t_results_queue* q);
int t_pool_results_queue_empty (t_results_queue* q);
int t_pool_results_queue_len (t_results_queue* q);
int t_pool_results_queue_sz (t_results_queue* q);
kh_m_i2i_t* kh_init_m_i2i ();
void kh_destroy_m_i2i (kh_m_i2i_t* h);
void kh_clear_m_i2i (kh_m_i2i_t* h);
khint_t kh_get_m_i2i (const(kh_m_i2i_t)* h, khint32_t key);
int kh_resize_m_i2i (kh_m_i2i_t* h, khint_t new_n_buckets);
khint_t kh_put_m_i2i (kh_m_i2i_t* h, khint32_t key, int* ret);
void kh_del_m_i2i (kh_m_i2i_t* h, khint_t x);
kh_s_i2i_t* kh_init_s_i2i ();
void kh_destroy_s_i2i (kh_s_i2i_t* h);
void kh_clear_s_i2i (kh_s_i2i_t* h);
khint_t kh_get_s_i2i (const(kh_s_i2i_t)* h, khint32_t key);
int kh_resize_s_i2i (kh_s_i2i_t* h, khint_t new_n_buckets);
khint_t kh_put_s_i2i (kh_s_i2i_t* h, khint32_t key, int* ret);
void kh_del_s_i2i (kh_s_i2i_t* h, khint_t x);
kh_map_t* kh_init_map ();
void kh_destroy_map (kh_map_t* h);
void kh_clear_map (kh_map_t* h);
khint_t kh_get_map (const(kh_map_t)* h, kh_cstr_t key);
int kh_resize_map (kh_map_t* h, khint_t new_n_buckets);
khint_t kh_put_map (kh_map_t* h, kh_cstr_t key, int* ret);
void kh_del_map (kh_map_t* h, khint_t x);
kh_refs_t* kh_init_refs ();
void kh_destroy_refs (kh_refs_t* h);
void kh_clear_refs (kh_refs_t* h);
khint_t kh_get_refs (const(kh_refs_t)* h, kh_cstr_t key);
int kh_resize_refs (kh_refs_t* h, khint_t new_n_buckets);
khint_t kh_put_refs (kh_refs_t* h, kh_cstr_t key, int* ret);
void kh_del_refs (kh_refs_t* h, khint_t x);
int is_directory (char* fn);
int is_file (char* fn);
int file_size (char* fn);
int itf8_decode (cram_fd* fd, int* val);
int itf8_put_blk (cram_block* blk, int val);
cram_block* cram_new_block (cram_content_type content_type, int content_id);
cram_block* cram_read_block (cram_fd* fd);
int cram_write_block (cram_fd* fd, cram_block* b);
void cram_free_block (cram_block* b);
char* zlib_mem_inflate (char* cdata, size_t csize, size_t* size);
int cram_uncompress_block (cram_block* b);
int cram_compress_block (cram_fd* fd, cram_block* b, cram_metrics* metrics, int method, int level);
cram_metrics* cram_new_metrics ();
char* cram_block_method2str (cram_block_method m);
char* cram_content_type2str (cram_content_type t);
ubyte* append_uint (ubyte* cp, int i);
int cram_load_reference (cram_fd* fd, const(char)* fn);
int refs2id (refs_t* r, SAM_hdr* bfd);
void refs_free (refs_t* r);
char* cram_get_ref (cram_fd* fd, int id, int start, int end);
void cram_ref_incr (refs_t* r, int id);
void cram_ref_decr (refs_t* r, int id);
cram_container* cram_new_container (int nrec, int nslice);
void cram_free_container (cram_container* c);
cram_container* cram_read_container (cram_fd* fd);
int cram_write_container (cram_fd* fd, cram_container* h);
int cram_flush_container (cram_fd* fd, cram_container* c);
int cram_flush_container_mt (cram_fd* fd, cram_container* c);
cram_block_compression_hdr* cram_new_compression_header ();
void cram_free_compression_header (cram_block_compression_hdr* hdr);
void cram_free_slice_header (cram_block_slice_hdr* hdr);
void cram_free_slice (cram_slice* s);
cram_slice* cram_new_slice (cram_content_type type, int nrecs);
cram_slice* cram_read_slice (cram_fd* fd);
cram_file_def* cram_read_file_def (cram_fd* fd);
int cram_write_file_def (cram_fd* fd, cram_file_def* def);
void cram_free_file_def (cram_file_def* def);
SAM_hdr* cram_read_SAM_hdr (cram_fd* fd);
int cram_write_SAM_hdr (cram_fd* fd, SAM_hdr* hdr);
cram_fd* cram_open (const(char)* filename, const(char)* mode);
cram_fd* cram_dopen (hFILE* fp, const(char)* filename, const(char)* mode);
int cram_close (cram_fd* fd);
int cram_seek (cram_fd* fd, off_t offset, int whence);
int cram_flush (cram_fd* fd);
int cram_eof (cram_fd* fd);
int cram_set_option (cram_fd* fd, cram_option opt, ...);
int cram_set_voption (cram_fd* fd, cram_option opt, va_list args);
int cram_set_header (cram_fd* fd, SAM_hdr* hdr);
int cram_put_bam_seq (cram_fd* fd, bam_seq_t* b);
cram_block* cram_encode_compression_header (cram_fd* fd, cram_container* c, cram_block_compression_hdr* h);
cram_block* cram_encode_slice_header (cram_fd* fd, cram_slice* s);
int cram_encode_container (cram_fd* fd, cram_container* c);
cram_record* cram_get_seq (cram_fd* fd);
int cram_get_bam_seq (cram_fd* fd, bam_seq_t** bam);
cram_block_compression_hdr* cram_decode_compression_header (cram_fd* fd, cram_block* b);
cram_block_slice_hdr* cram_decode_slice_header (cram_fd* fd, cram_block* b);
int cram_decode_slice (cram_fd* fd, cram_container* c, cram_slice* s, SAM_hdr* hdr);
cram_stats* cram_stats_create ();
void cram_stats_add (cram_stats* st, int val);
void cram_stats_del (cram_stats* st, int val);
void cram_stats_dump (cram_stats* st);
void cram_stats_free (cram_stats* st);
cram_encoding cram_stats_encoding (cram_fd* fd, cram_stats* st);
char* cram_encoding2str (cram_encoding t);
cram_codec* cram_decoder_init (cram_encoding codec, char* data, int size, cram_external_type option, int version_);
cram_codec* cram_encoder_init (cram_encoding codec, cram_stats* st, cram_external_type option, void* dat, int version_);
int cram_codec_to_id (cram_codec* c, int* id2);
int cram_index_load (cram_fd* fd, const(char)* fn);
void cram_index_free (cram_fd* fd);
cram_index* cram_index_query (cram_fd* fd, int refid, int pos, cram_index* frm);
int cram_seek_to_refpos (cram_fd* fd, cram_range* r);
void cram_index_free (cram_fd* fd);
int cram_index_build (cram_fd* fd, const(char)* fn_base);
int cram_to_bam(SAM_hdr*, cram_fd*, cram_slice*, cram_record*, int, bam_seq_t**);
