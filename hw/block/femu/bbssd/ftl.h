#ifndef __FEMU_FTL_H
#define __FEMU_FTL_H

#define FEMU_DEBUG_FTL

#include "../nvme.h"

#define INVALID_PPA     (~(0ULL))
#define INVALID_LPN     (~(0ULL))
#define UNMAPPED_PPA    (~(0ULL))


// DZ Start
// Define I/O operations
#define INITIAL_OP    0 // Supposed to be flush, but not used by SNAP. Consider as record not used before.
#define WRITE_OP      1 // Write
#define WRITE_ONCE_OP 2 // Only wrote once, will segue into WRITE_OP on next write
#define DISCARD_OP    3 // Discard

// Define decay factor
#define DECAY 0.9

// Time Precision in Bits
#define TIME_PREC_BITS    (15)
// ms to s
#define MS_PER_S 1000


#define USE_FREE  0
#define USE_INUSE 1
#define USE_FULL  2
// DZ End

enum {
    NAND_READ =  0,
    NAND_WRITE = 1,
    NAND_ERASE = 2,

    NAND_READ_LATENCY = 40000,
    NAND_PROG_LATENCY = 200000,
    NAND_ERASE_LATENCY = 2000000,
    // DZ Start
    // Define fast programming latency
    NAND_FAST_PROG_LATENCY = 150000,
    // DZ End
};

enum {
    USER_IO = 0,
    GC_IO = 1,
};

enum {
    SEC_FREE = 0,
    SEC_INVALID = 1,
    SEC_VALID = 2,

    PG_FREE = 0,
    PG_INVALID = 1,
    PG_VALID = 2
};

enum {
    FEMU_ENABLE_GC_DELAY = 1,
    FEMU_DISABLE_GC_DELAY = 2,

    FEMU_ENABLE_DELAY_EMU = 3,
    FEMU_DISABLE_DELAY_EMU = 4,

    FEMU_RESET_ACCT = 5,
    FEMU_ENABLE_LOG = 6,
    FEMU_DISABLE_LOG = 7,
    FEMU_SHOW_WA_INFO = 8,
    FEMU_LOG_FLUSH = 9,
};


#define BLK_BITS    (16)
#define PG_BITS     (16)
#define SEC_BITS    (8)
#define PL_BITS     (8)
#define LUN_BITS    (8)
#define CH_BITS     (7)

/* describe a physical page addr */
struct ppa {
    union {
        struct {
            uint64_t blk : BLK_BITS;
            uint64_t pg  : PG_BITS;
            uint64_t sec : SEC_BITS;
            uint64_t pl  : PL_BITS;
            uint64_t lun : LUN_BITS;
            uint64_t ch  : CH_BITS;
            uint64_t rsv : 1;
        } g;

        uint64_t ppa;
    };
};

typedef int nand_sec_status_t;

struct nand_page {
    nand_sec_status_t *sec;
    int nsecs;
    int status;
};

struct nand_block {
    struct nand_page *pg;
    int npgs;
    int ipc; /* invalid page count */
    int vpc; /* valid page count */
    int erase_cnt;
    int wp; /* current write pointer */
    int status;
};

struct nand_plane {
    struct nand_block *blk;
    int nblks;
};

struct nand_lun {
    struct nand_plane *pl;
    int npls;
    uint64_t next_lun_avail_time;
    bool busy;
    uint64_t gc_endtime;
};

struct ssd_channel {
    struct nand_lun *lun;
    int nluns;
    uint64_t next_ch_avail_time;
    bool busy;
    uint64_t gc_endtime;
};

struct ssdparams {
    int secsz;        /* sector size in bytes */
    int secs_per_pg;  /* # of sectors per page */
    int pgs_per_blk;  /* # of NAND pages per block */
    int blks_per_pl;  /* # of blocks per plane */
    int pls_per_lun;  /* # of planes per LUN (Die) */
    int luns_per_ch;  /* # of LUNs per channel */
    int nchs;         /* # of channels in the SSD */

    int pg_rd_lat;    /* NAND page read latency in nanoseconds */
    int pg_wr_lat;    /* NAND page program latency in nanoseconds */
    int blk_er_lat;   /* NAND block erase latency in nanoseconds */
    int ch_xfer_lat;  /* channel transfer latency for one page in nanoseconds
                       * this defines the channel bandwith
                       */

    double gc_thres_pcent;
    int gc_thres_lines;
    int gc_thres_blocks;
    double gc_thres_pcent_high;
    int gc_thres_lines_high;
    int gc_thres_blocks_high;
    bool enable_gc_delay;

    /* below are all calculated values */
    int secs_per_blk; /* # of sectors per block */
    int secs_per_pl;  /* # of sectors per plane */
    int secs_per_lun; /* # of sectors per LUN */
    int secs_per_ch;  /* # of sectors per channel */
    int tt_secs;      /* # of sectors in the SSD */

    int pgs_per_pl;   /* # of pages per plane */
    int pgs_per_lun;  /* # of pages per LUN (Die) */
    int pgs_per_ch;   /* # of pages per channel */
    int tt_pgs;       /* total # of pages in the SSD */

    int min_blks_per_lun; /* min # of blocks per LUN */
    int blks_per_lun; /* # of blocks per LUN, supports using several of all channels */
    int blks_per_ch;  /* # of blocks per channel */
    int tt_blks;      /* total # of blocks in the SSD */

    int secs_per_line;
    int pgs_per_line;
    int blks_per_line;
    int tt_lines;

    int pls_per_ch;   /* # of planes per channel */
    int tt_pls;       /* total # of planes in the SSD */

    int tt_luns;      /* total # of LUNs in the SSD */

    int tt_chunks;    /* total # of chunks in the SSD */

    int pages_per_superblock;

    bool death_time_prediction; /* Death Time Prediction Enabled */
    bool enable_stream;
    bool enable_stream_redirect;
    uint32_t channel_split_exp;
    uint32_t default_channels_per_line;
    uint32_t default_luns_per_channel;
    uint32_t init_blk_per_plane;
    int channel_regions;
    int lun_regions;
    int blk_regions;
    uint8_t msl;
    uint64_t epoch;   /* Last updated age */
    uint64_t start_time; /* When the system started */
    uint64_t max_age;
    uint32_t access_interval_precision;

    int real_num_streams;
    int gc_stream_id;
    bool enable_hetero_sbsize;
    bool enable_partial_gc;
    int gc_start_channel;
    int gc_start_lun;

    int max_allow_gc_lines;
};

typedef struct line {
    QTAILQ_ENTRY(line) **full_entry; /* in either {free,victim,full} list */
    QTAILQ_ENTRY(line) **victim_entry; /* in either {free,victim,full} list */
    int id;  /* line id, the same as corresponding block id */
    int ipc; /* invalid page count in this line */
    int vpc; /* valid page count in this line */
    int pgs_per_line;  /* pages in this line */
    int start_channel; /* The start channel # for this line to use */
    int start_lun; /* The start channel # for this line to use */
    int *channel_list; /* The list of channels IDs to use */
    int **lun_list;     /* The list of LUNs to use */
    int **block_list;   /* The list of block IDs to use */
    int total_channels; /* Number of channels this line can use */
    int *total_luns;    /* Number of LUNs per channel used */
    
    
    bool **inserted_to_full_queue;
    bool **inserted_to_victim_queue;
    bool valid;
    int stream;
    double earliest_dt;
    double latest_dt;
    double close_time;
    double expected_h;
    int use;
    int curchannel;
    int curlun;
} line;

/* wp: record next write addr */
struct write_pointer {
    struct line *curline;
    int ch;
    int lun;
    int pg;
    int blk;
    int pl;
};

struct stream_info {
    // The following entries are defined for outgoing block
    uint64_t earliest_death_time;
    uint64_t latest_death_time;
    uint64_t block_open_time;
    bool sender;
    bool receiver;


    // The following entries are defined for the incoming stream
    bool full_before;
    double avg_incoming_interval;
    double avg_temp_incoming_interval;
    uint64_t stream_counter_start_time;
    uint64_t next_avail_time;
    double avg_pages_per_request;
};

struct line_resource_mgmt{
    /* free line list, we only need to maintain a list of blk numbers */
    QTAILQ_HEAD(free_line_list, line) free_line_list;
    QTAILQ_HEAD(victim_line_list, line) victim_line_list;
    QTAILQ_HEAD(full_line_list, line) full_line_list;
    int free_line_cnt;
    int victim_line_cnt;
    int full_line_cnt;
};

struct line_mgmt {
    struct line *lines;
    struct line_resource_mgmt** line_resource;
    int tt_lines;
    struct line ****channel_lines;
};

struct nand_cmd {
    int type;
    int cmd;
    int64_t stime; /* Coperd: request arrival time */
};

// DZ Start
/*
struct death_time_track{
    // Define structure for death time analysis
    // Previously calculated average
    uint64_t death_time_avg;
    #ifdef FEMU_DEBUG_FTL
    uint64_t prev_death_time_prediction;
    #endif
    // Previous access timestamp
    uint64_t age;
    // Previous access I/O operation
    // Can be READ_OP, WRITE_OP, FLUSH_OP or DISCARD_OP
    int last_access_op;
    bool valid;
};
*/

struct wa_info {
    uint64_t pages_from_host;
    uint64_t pages_from_gc;
};

struct death_time_track {
    uint64_t death_time_avg    : TIME_PREC_BITS;
    uint64_t age               : TIME_PREC_BITS;
    uint64_t last_access_op    : 2;
    #ifdef FEMU_DEBUG_FTL
    uint64_t prev_death_time_prediction : TIME_PREC_BITS;
    #endif
};

struct seq_write_info {
    uint size;
    uint cur;
    uint64_t *list;
    uint inited;
};

struct block_num{
    int block_num;
    QTAILQ_ENTRY(block_num) entry;
};

struct block_mgmt {
    int free_blocks_cnt;
    QTAILQ_HEAD(free_block_list, block_num) free_block_list;
};

struct lun_mgmt{
    int next_avail_lun;
    struct block_mgmt* lun;
};

struct channel_mgmt{
    int next_avail_channel;
    uint64_t next_line_id;
    struct lun_mgmt* channel;
};
// DZ End

struct ssd {
    char *ssdname;
    struct ssdparams sp;
    struct ssd_channel *ch;
    struct ppa *maptbl; /* page level mapping table */
    // DZ Start
    // We have a structure for death time analysis, which splits LBA into chunks
    // Number of pages/chunk is defined as pages_per_chunk in FEMU start script
    struct death_time_track *death_time_list; /* page level mapping table */
    struct stream_info *stream_info;
    // DZ End
    uint64_t *rmap;     /* reverse mapptbl, assume it's stored in OOB */
    struct write_pointer *wp; // We need multiple pointers for multi-stream SSD.
    struct line_mgmt lm;

    /* lockless ring for communication with NVMe IO thread */
    struct rte_ring **to_ftl;
    struct rte_ring **to_poller;
    bool *dataplane_started_ptr;
    QemuThread ftl_thread;
    uint64_t pages_from_host;
    uint64_t pages_from_gc;
    struct seq_write_info *seq_info;
    struct channel_mgmt channel_mgmt;
};

void ssd_init(FemuCtrl *n);

uint64_t get_uptime(struct ssd *ssd);

#ifdef FEMU_DEBUG_FTL
#define ftl_debug(fmt, ...) \
    do { printf("[FEMU] FTL-Dbg: " fmt, ## __VA_ARGS__); } while (0)
#else
#define ftl_debug(fmt, ...) \
    do { } while (0)
#endif

#define ftl_err(fmt, ...) \
    do { fprintf(stderr, "[FEMU] FTL-Err: " fmt, ## __VA_ARGS__); } while (0)

#define ftl_log(fmt, ...) \
    do { printf("[FEMU] FTL-Log: " fmt, ## __VA_ARGS__); } while (0)


/* FEMU assert() */
#ifdef FEMU_DEBUG_FTL
#define ftl_assert(expression) \
    do {fflush(NULL); assert(expression);} while (0)
#else
#define ftl_assert(expression)
#endif

#endif
