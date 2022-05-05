#include "./ftl.h"
#include <math.h>
#ifdef FEMU_DEBUG_FTL
FILE * femu_log_file;
#define write_log(fmt, ...) \
    do { if (femu_log_file) fprintf(femu_log_file, fmt, ## __VA_ARGS__);} while (0)
#else
#define write_log(fmt, ...) \
    do { } while (0)
#endif

static void *ftl_thread(void *arg);

static inline bool should_gc(struct ssd *ssd)
{
    return (ssd->lm.free_line_cnt <= ssd->sp.gc_thres_lines);
}

static inline bool should_gc_high(struct ssd *ssd)
{
    return (ssd->lm.free_line_cnt <= ssd->sp.gc_thres_lines_high);
}

static inline struct ppa get_maptbl_ent(struct ssd *ssd, uint64_t lpn)
{
    return ssd->maptbl[lpn];
}

static inline void set_maptbl_ent(struct ssd *ssd, uint64_t lpn, struct ppa *ppa)
{
    ftl_assert(lpn < ssd->sp.tt_pgs);
    ssd->maptbl[lpn] = *ppa;
}

// DZ Start
static uint64_t lpn_to_chunk(uint64_t lpn, uint32_t pages_per_chunk){
    // TODO: add hashing support
    return lpn / pages_per_chunk;
}

static uint64_t get_chunk_list(FemuCtrl *n, struct ssd *ssd, uint64_t start_lpn, uint64_t end_lpn, uint64_t* chunks){
    uint64_t i;
    //write_log("Calculating chunks from LPN %"PRIu64" to LPN %"PRIu64"\n", start_lpn, end_lpn);
    if (n->enable_hashing){
        for (i = start_lpn; i <= end_lpn; i++){
            chunks[i - start_lpn] = lpn_to_chunk(i, n->pages_per_chunk);
            //write_log("LPN %"PRIu64" has chunk num: %"PRIu64"\n", i, chunks[i - start_lpn]);
        }
        return (end_lpn - start_lpn) + 1;
    }else{
        uint64_t start_chunk = lpn_to_chunk(start_lpn, n->pages_per_chunk);
        //write_log("Start LPN %"PRIu64" has chunk num: %"PRIu64"\n", start_lpn, start_chunk);
        uint64_t end_chunk = lpn_to_chunk(end_lpn, n->pages_per_chunk);
        //write_log("End LPN %"PRIu64" has chunk num: %"PRIu64"\n", end_lpn, end_chunk);
        for (i = start_chunk; i <= end_chunk; i++){
            chunks[i - start_chunk] = i;
            //write_log("Chunk i=%"PRIu64" has chunk num: %"PRIu64"\n", i - start_chunk, chunks[i - start_chunk]);
        }
        return (end_chunk - start_chunk) + 1;
    }
}

static uint64_t get_passed_epoch_since_update(struct ssd *ssd){
    uint64_t now_real_time = qemu_clock_get_ms(QEMU_CLOCK_REALTIME);
    uint64_t passed_epoch = (now_real_time - ssd->sp.epoch) / MS_PER_S / ssd->sp.access_interval_precision;
    return passed_epoch;
}

static uint64_t get_passed_epoch_since_start(struct ssd *ssd){
    uint64_t now_real_time = qemu_clock_get_ms(QEMU_CLOCK_REALTIME);
    uint64_t system_epoch = (now_real_time - ssd->sp.start_time) / MS_PER_S / ssd->sp.access_interval_precision;
    return system_epoch;
}

uint64_t get_uptime(struct ssd *ssd){
    uint64_t now_real_time = qemu_clock_get_ms(QEMU_CLOCK_REALTIME);
    uint64_t uptime = (now_real_time - ssd->sp.start_time) / MS_PER_S;
    return uptime;
}

static void set_latest_access_time(FemuCtrl *n, struct ssd *ssd, uint64_t start_lpn, uint64_t end_lpn, int op)
{
    uint64_t passed_epoch = get_passed_epoch_since_update(ssd);
    if (passed_epoch > 0){
        // Update the last updated timestamp
        ssd->sp.epoch += ssd->sp.access_interval_precision * MS_PER_S * passed_epoch;
        // Update ages of all pages
        for (int i = 0; i < ssd->sp.tt_chunks; i++){
            if ((ssd->death_time_list[i].age + passed_epoch > ssd->sp.max_age)){
                ssd->death_time_list[i].age = ssd->sp.max_age;
            }else{
                ssd->death_time_list[i].age += passed_epoch;
            }
        }
    }

    int      prev_op;
    uint64_t prev_avg;
    uint64_t prev_age;
    uint64_t i;
    uint64_t chunk;

    // If hashing is not enabled, we can just skip to next chunk's page
    // So that the same chunk will not be reupdated
    uint64_t *chunks = g_malloc0(sizeof(uint64_t) * ((end_lpn - start_lpn) + 1));
    uint64_t num_chunks = get_chunk_list(n, ssd, start_lpn, end_lpn, chunks);
    
    for (i = 0; i < num_chunks; i++){
        chunk = chunks[i];
        if (ssd->death_time_list[chunk].last_access_op != INITIAL_OP){
            prev_op = ssd->death_time_list[chunk].last_access_op;
            prev_avg = ssd->death_time_list[chunk].death_time_avg;
            prev_age = ssd->death_time_list[chunk].age;
            
            // Only consider W->W and W->D as death. Update death time avg in this case
            if (prev_op == WRITE_OP || prev_op == WRITE_ONCE_OP){
                if (prev_avg > 0){
                    ssd->death_time_list[chunk].death_time_avg = prev_avg * (1 - DECAY) + prev_age * DECAY;
                    //write_log("Old prediction: %"PRIu64", age: %"PRIu64", new prediction: %d, Chunk: %"PRIu64"\n", prev_avg, prev_age, ssd->death_time_list[chunk].death_time_avg, chunk);
                }else{
                    if (likely(prev_op == WRITE_OP)){
                        ssd->death_time_list[chunk].death_time_avg = prev_avg * (1 - DECAY) + prev_age * DECAY;
                        //write_log("Old prediction: %"PRIu64", age: %"PRIu64", new prediction: %d, Chunk: %"PRIu64"\n", prev_avg, prev_age, ssd->death_time_list[chunk].death_time_avg, chunk);
                    }else{
                        ftl_assert(prev_avg == 0);
                        //write_log("First death, age: %"PRIu64", new prediction: %"PRIu64", Chunk: %"PRIu64"\n", prev_age, prev_age, chunk);
                        ssd->death_time_list[chunk].death_time_avg = prev_age;
                    }
                }
            }else{
                #ifdef FEMU_DEBUG_FTL
                ////write_log("Discard for chunk: %"PRIu64"\n", chunk);
                #endif
            }

            #ifdef FEMU_DEBUG_FTL
            ssd->death_time_list[chunk].prev_death_time_prediction = prev_avg;
            #endif
            ssd->death_time_list[chunk].last_access_op = op;
            // Clear age, since this block is now dead
            ssd->death_time_list[chunk].age = 0;
        }else{
            if (op == WRITE_OP){
                chunk = chunks[i];
                ssd->death_time_list[chunk].death_time_avg = 0;
                #ifdef FEMU_DEBUG_FTL
                ssd->death_time_list[chunk].prev_death_time_prediction = 0;
                //write_log("First time write at chunk: %"PRIu64"\n", chunk);
                #endif
                ssd->death_time_list[chunk].last_access_op = WRITE_ONCE_OP;
                ssd->death_time_list[chunk].age = 0;
            }else{
                //write_log("Discarding never-used chunk: %"PRIu64"\n", chunk);
            }
        }
    }
}
// DZ End

static uint64_t ppa2pgidx(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    uint64_t pgidx;

    pgidx = ppa->g.ch  * spp->pgs_per_ch  + \
            ppa->g.lun * spp->pgs_per_lun + \
            ppa->g.pl  * spp->pgs_per_pl  + \
            ppa->g.blk * spp->pgs_per_blk + \
            ppa->g.pg;

    ftl_assert(pgidx < spp->tt_pgs);

    return pgidx;
}

static inline uint64_t get_rmap_ent(struct ssd *ssd, struct ppa *ppa)
{
    uint64_t pgidx = ppa2pgidx(ssd, ppa);

    return ssd->rmap[pgidx];
}

/* set rmap[page_no(ppa)] -> lpn */
static inline void set_rmap_ent(struct ssd *ssd, uint64_t lpn, struct ppa *ppa)
{
    uint64_t pgidx = ppa2pgidx(ssd, ppa);

    ssd->rmap[pgidx] = lpn;
}

static inline int victim_line_cmp_pri(pqueue_pri_t next, pqueue_pri_t curr)
{
    return (next > curr);
}

static inline pqueue_pri_t victim_line_get_pri(void *a)
{
    return ((struct line *)a)->vpc;
}

static inline void victim_line_set_pri(void *a, pqueue_pri_t pri)
{
    ((struct line *)a)->vpc = pri;
}

static inline size_t victim_line_get_pos(void *a)
{
    return ((struct line *)a)->pos;
}

static inline void victim_line_set_pos(void *a, size_t pos)
{
    ((struct line *)a)->pos = pos;
}

static void ssd_init_lines(struct ssd *ssd)
{
    struct ssdparams *spp = &ssd->sp;
    struct line_mgmt *lm = &ssd->lm;
    struct line *line;

    lm->tt_lines = spp->blks_per_pl;
    ftl_assert(lm->tt_lines == spp->tt_lines);
    lm->lines = g_malloc0(sizeof(struct line) * lm->tt_lines);

    QTAILQ_INIT(&lm->free_line_list);
    lm->victim_line_pq = pqueue_init(spp->tt_lines, victim_line_cmp_pri,
            victim_line_get_pri, victim_line_set_pri,
            victim_line_get_pos, victim_line_set_pos);
    QTAILQ_INIT(&lm->full_line_list);

    lm->free_line_cnt = 0;
    for (int i = 0; i < lm->tt_lines; i++) {
        line = &lm->lines[i];
        line->id = i;
        line->ipc = 0;
        line->vpc = 0;
        line->pos = 0;
        /* initialize all the lines as free lines */
        QTAILQ_INSERT_TAIL(&lm->free_line_list, line, entry);
        lm->free_line_cnt++;
    }

    ftl_assert(lm->free_line_cnt == lm->tt_lines);
    lm->victim_line_cnt = 0;
    lm->full_line_cnt = 0;
}

static void ssd_init_write_pointer(struct ssd *ssd, uint8_t streams)
{
    // n streams -> n+1 write pointers since we need stream 0 as default stream.
    ssd->wp = g_malloc0(sizeof(struct write_pointer) * (streams + 1));
    ssd->stream_info = g_malloc0(sizeof(struct stream_info) * (streams + 1));
    struct write_pointer *wpp = NULL;
    struct stream_info *si = NULL;
    struct line_mgmt *lm = &ssd->lm;
    struct line *curline = NULL;
    for (int i = 0; i <= streams; i++){
        curline = QTAILQ_FIRST(&lm->free_line_list);
        QTAILQ_REMOVE(&lm->free_line_list, curline, entry);
        lm->free_line_cnt--;
        wpp = &ssd->wp[i];
        si = &ssd->stream_info[i];
        /* wpp->curline is always our next-to-write super-block */
        wpp->curline = curline;
        wpp->ch = 0;
        wpp->lun = 0;
        wpp->pg = 0;
        wpp->blk = i;
        wpp->pl = 0;
        si->earliest_death_time = 0;
        si->latest_death_time = 0;
        si->avg_incoming_interval = 0;
        si->block_start_time = get_uptime(ssd);
        si->fulled_before = false;
        si->sender = false;
        si->receiver = false;
        si->page_counter = 0;
    }
}

static inline void check_addr(int a, int max)
{
    ftl_assert(a >= 0 && a < max);
}

static struct line *get_next_free_line(struct ssd *ssd)
{
    struct line_mgmt *lm = &ssd->lm;
    struct line *curline = NULL;

    curline = QTAILQ_FIRST(&lm->free_line_list);
    if (!curline) {
        ftl_err("No free lines left in [%s] !!!!\n", ssd->ssdname);
        return NULL;
    }

    QTAILQ_REMOVE(&lm->free_line_list, curline, entry);
    lm->free_line_cnt--;
    return curline;
}

static void ssd_advance_write_pointer(struct ssd *ssd, uint8_t stream)
{
    struct ssdparams *spp = &ssd->sp;
    struct write_pointer *wpp = &ssd->wp[stream];
    struct stream_info *si = &ssd->stream_info[stream];
    si->page_counter++;
    uint64_t uptime = get_uptime(ssd);
    struct line_mgmt *lm = &ssd->lm;
    check_addr(wpp->ch, spp->nchs);
    wpp->ch++;
    if (wpp->ch == spp->nchs) {
        wpp->ch = 0;
        check_addr(wpp->lun, spp->luns_per_ch);
        wpp->lun++;
        /* in this case, we should go to next lun */
        if (wpp->lun == spp->luns_per_ch) {
            wpp->lun = 0;
            /* go to next page in the block */
            check_addr(wpp->pg, spp->pgs_per_blk);
            wpp->pg++;
            if (wpp->pg == spp->pgs_per_blk) {
                wpp->pg = 0;
                /* move current line to {victim,full} line list */
                if (wpp->curline->vpc == spp->pgs_per_line) {
                    /* all pgs are still valid, move to full line list */
                    ftl_assert(wpp->curline->ipc == 0);
                    QTAILQ_INSERT_TAIL(&lm->full_line_list, wpp->curline, entry);
                    lm->full_line_cnt++;
                } else {
                    ftl_assert(wpp->curline->vpc >= 0 && wpp->curline->vpc < spp->pgs_per_line);
                    /* there must be some invalid pages in this line */
                    ftl_assert(wpp->curline->ipc > 0);
                    pqueue_insert(lm->victim_line_pq, wpp->curline);
                    lm->victim_line_cnt++;
                }
                /* current line is used up, pick another empty line */

                // DZ Start

                if (spp->enable_stream_redirect && spp->death_time_prediction){
                    //write_log("Stream %d block is full.\n", stream);
                    if (si->fulled_before){
                        si->avg_incoming_interval = si->avg_incoming_interval * (1 - DECAY) + ((uptime - si->block_start_time) / spp->pages_per_superblock) * DECAY;
                    }else{
                        //write_log("Stream %d first full\n", stream);
                        si->avg_incoming_interval = (uptime - si->block_start_time / spp->pages_per_superblock);
                        //write_log("Stream %d new avg full interval: %"PRIu64"\n", stream, si->avg_incoming_interval);
                        si->fulled_before = true;
                    }
                    
                }else{
                    // Prediction not used. Just skip rotation.
                }
                write_log("+++++\n\n");
                write_log("Stream %d curblock is full: \n", stream);
                write_log("Avg incoming interval: %.6fs\n", si->avg_incoming_interval);
                write_log("Block open time: %"PRIu64"s\n", si->block_start_time);
                write_log("Block close time: %"PRIu64"s\n", uptime);
                write_log("Block earliest DT: %"PRIu64" units\n", si->earliest_death_time);
                write_log("Block latest DT: %"PRIu64" units\n", si->latest_death_time);
                write_log("Transient time: %"PRIu64" units\n", si->latest_death_time - (get_passed_epoch_since_start(ssd) > si->earliest_death_time ? get_passed_epoch_since_start(ssd) : si->earliest_death_time));
                write_log("-----\n\n");
                si->sender = false;
                si->receiver = false;
                si->block_start_time = uptime;
                si->earliest_death_time = 0;
                si->latest_death_time = 0;
                ftl_assert(si->page_counter == spp->pages_per_superblock);
                si->page_counter = 0;

                check_addr(wpp->blk, spp->blks_per_pl);
                wpp->curline = NULL;
                wpp->curline = get_next_free_line(ssd);
                if (!wpp->curline) {
                    /* TODO */
                    abort();
                }
                wpp->blk = wpp->curline->id;
                //femu_log("New superblock: %d caused by stream: %d\n", wpp->blk, stream);
                check_addr(wpp->blk, spp->blks_per_pl);
                /* make sure we are starting from page 0 in the super block */
                ftl_assert(wpp->pg == 0);
                ftl_assert(wpp->lun == 0);
                ftl_assert(wpp->ch == 0);
                /* TODO: assume # of pl_per_lun is 1, fix later */
                ftl_assert(wpp->pl == 0);
            }
        }
    }
}

static struct ppa get_new_page(struct ssd *ssd, uint8_t stream)
{
    struct write_pointer *wpp = &ssd->wp[stream];
    struct ppa ppa;
    ppa.ppa = 0;
    ppa.g.ch = wpp->ch;
    ppa.g.lun = wpp->lun;
    ppa.g.pg = wpp->pg;
    ppa.g.blk = wpp->blk;
    ppa.g.pl = wpp->pl;
    //printf("Stream: %d, channel: %d, lun: %d, pg: %d, blk: %d, Plane: %d\n", stream, ppa.g.ch, ppa.g.lun, ppa.g.pg, ppa.g.blk, ppa.g.pl);
    //fflush(NULL);
    ftl_assert(ppa.g.pl == 0);

    return ppa;
}

static void check_params(struct ssdparams *spp)
{
    /*
     * we are using a general write pointer increment method now, no need to
     * force luns_per_ch and nchs to be power of 2
     */

    //ftl_assert(is_power_of_2(spp->luns_per_ch));
    //ftl_assert(is_power_of_2(spp->nchs));
}

static void ssd_init_params(struct ssdparams *spp)
{
    spp->secsz = 512;
    spp->secs_per_pg = 8;
    spp->pgs_per_blk = 256;
    spp->blks_per_pl = 320; /* 20GB */
    spp->pls_per_lun = 1;
    spp->luns_per_ch = 8;
    spp->nchs = 8;

    spp->pg_rd_lat = NAND_READ_LATENCY;
    spp->pg_wr_lat = NAND_PROG_LATENCY;
    spp->blk_er_lat = NAND_ERASE_LATENCY;
    spp->ch_xfer_lat = 0;

    /* calculated values */
    spp->secs_per_blk = spp->secs_per_pg * spp->pgs_per_blk;
    spp->secs_per_pl = spp->secs_per_blk * spp->blks_per_pl;
    spp->secs_per_lun = spp->secs_per_pl * spp->pls_per_lun;
    spp->secs_per_ch = spp->secs_per_lun * spp->luns_per_ch;
    spp->tt_secs = spp->secs_per_ch * spp->nchs;

    spp->pgs_per_pl = spp->pgs_per_blk * spp->blks_per_pl;
    spp->pgs_per_lun = spp->pgs_per_pl * spp->pls_per_lun;
    spp->pgs_per_ch = spp->pgs_per_lun * spp->luns_per_ch;
    spp->tt_pgs = spp->pgs_per_ch * spp->nchs;

    spp->blks_per_lun = spp->blks_per_pl * spp->pls_per_lun;
    spp->blks_per_ch = spp->blks_per_lun * spp->luns_per_ch;
    spp->tt_blks = spp->blks_per_ch * spp->nchs;

    spp->pls_per_ch =  spp->pls_per_lun * spp->luns_per_ch;
    spp->tt_pls = spp->pls_per_ch * spp->nchs;

    spp->tt_luns = spp->luns_per_ch * spp->nchs;

    spp->pages_per_superblock = spp->pgs_per_blk * spp->luns_per_ch * spp->nchs;

    /* line is special, put it at the end */
    spp->blks_per_line = spp->tt_luns; /* TODO: to fix under multiplanes */
    spp->pgs_per_line = spp->blks_per_line * spp->pgs_per_blk;
    spp->secs_per_line = spp->pgs_per_line * spp->secs_per_pg;
    spp->tt_lines = spp->blks_per_lun; /* TODO: to fix under multiplanes */

    spp->gc_thres_pcent = 0.75;
    spp->gc_thres_lines = (int)((1 - spp->gc_thres_pcent) * spp->tt_lines);
    spp->gc_thres_pcent_high = 0.95;
    spp->gc_thres_lines_high = (int)((1 - spp->gc_thres_pcent_high) * spp->tt_lines);
    spp->enable_gc_delay = true;


    check_params(spp);
}

static void ssd_init_nand_page(struct nand_page *pg, struct ssdparams *spp)
{
    pg->nsecs = spp->secs_per_pg;
    pg->sec = g_malloc0(sizeof(nand_sec_status_t) * pg->nsecs);
    for (int i = 0; i < pg->nsecs; i++) {
        pg->sec[i] = SEC_FREE;
    }
    pg->status = PG_FREE;
}

static void ssd_init_nand_blk(struct nand_block *blk, struct ssdparams *spp)
{
    blk->npgs = spp->pgs_per_blk;
    blk->pg = g_malloc0(sizeof(struct nand_page) * blk->npgs);
    for (int i = 0; i < blk->npgs; i++) {
        ssd_init_nand_page(&blk->pg[i], spp);
    }
    blk->ipc = 0;
    blk->vpc = 0;
    blk->erase_cnt = 0;
    blk->wp = 0;
}

static void ssd_init_nand_plane(struct nand_plane *pl, struct ssdparams *spp)
{
    pl->nblks = spp->blks_per_pl;
    pl->blk = g_malloc0(sizeof(struct nand_block) * pl->nblks);
    for (int i = 0; i < pl->nblks; i++) {
        ssd_init_nand_blk(&pl->blk[i], spp);
    }
}

static void ssd_init_nand_lun(struct nand_lun *lun, struct ssdparams *spp)
{
    lun->npls = spp->pls_per_lun;
    lun->pl = g_malloc0(sizeof(struct nand_plane) * lun->npls);
    for (int i = 0; i < lun->npls; i++) {
        ssd_init_nand_plane(&lun->pl[i], spp);
    }
    lun->next_lun_avail_time = 0;
    lun->busy = false;
}

static void ssd_init_ch(struct ssd_channel *ch, struct ssdparams *spp)
{
    ch->nluns = spp->luns_per_ch;
    ch->lun = g_malloc0(sizeof(struct nand_lun) * ch->nluns);
    for (int i = 0; i < ch->nluns; i++) {
        ssd_init_nand_lun(&ch->lun[i], spp);
    }
    ch->next_ch_avail_time = 0;
    ch->busy = 0;
}

static void ssd_init_maptbl(struct ssd *ssd)
{
    struct ssdparams *spp = &ssd->sp;

    ssd->maptbl = g_malloc0(sizeof(struct ppa) * spp->tt_pgs);
    for (int i = 0; i < spp->tt_pgs; i++) {
        ssd->maptbl[i].ppa = UNMAPPED_PPA;
    }
}

static void ssd_init_death_time(struct ssd *ssd, uint32_t pages_per_chunk)
{
    
    struct ssdparams *spp = &ssd->sp;

    spp->tt_chunks = (spp->tt_pgs - 1) / pages_per_chunk + 1;
    spp->epoch = qemu_clock_get_ms(QEMU_CLOCK_REALTIME);
    spp->start_time = spp->epoch;
    //write_log("Starting system... spp->epoch is now %"PRIu64"\n", spp->epoch);
    spp->max_age = (1 << TIME_PREC_BITS) - 1;

    ssd->death_time_list = g_malloc0(sizeof(struct death_time_track) * spp->tt_chunks);
    for (int i = 0; i < spp->tt_chunks; i++) {
        ssd->death_time_list[i].last_access_op = INITIAL_OP;
        ssd->death_time_list[i].death_time_avg = 0;
    }
}

static void ssd_init_rmap(struct ssd *ssd)
{
    struct ssdparams *spp = &ssd->sp;

    ssd->rmap = g_malloc0(sizeof(uint64_t) * spp->tt_pgs);
    for (int i = 0; i < spp->tt_pgs; i++) {
        ssd->rmap[i] = INVALID_LPN;
    }
}

void ssd_init(FemuCtrl *n)
{
    struct ssd *ssd = n->ssd;
    ssd->pages_from_host = 0;
    ssd->pages_from_gc = 0;
    
    struct ssdparams *spp = &ssd->sp;

    ftl_assert(ssd);

    ssd_init_params(spp);

    /* initialize ssd internal layout architecture */
    ssd->ch = g_malloc0(sizeof(struct ssd_channel) * spp->nchs);
    for (int i = 0; i < spp->nchs; i++) {
        ssd_init_ch(&ssd->ch[i], spp);
    }

    // Pass if we have enabled death time analysis
    spp->death_time_prediction = n->death_time_prediction;

    // Pass number of streams supported
    spp->enable_stream = n->enable_stream;
    spp->msl = n->msl;
    spp->enable_stream_redirect = n->enable_stream_redirect;
    spp->access_interval_precision = n->access_interval_precision;

    /* initialize maptbl */
    ssd_init_maptbl(ssd);

    /* initialize death_time */
    ssd_init_death_time(ssd, n->pages_per_chunk);

    /* initialize rmap */
    ssd_init_rmap(ssd);

    /* initialize all the lines */
    ssd_init_lines(ssd);

    /* initialize write pointer, this is how we allocate new pages for writes */
    ssd_init_write_pointer(ssd, n->msl);

    qemu_thread_create(&ssd->ftl_thread, "FEMU-FTL-Thread", ftl_thread, n,
                       QEMU_THREAD_JOINABLE);
}

static inline bool valid_ppa(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    int ch = ppa->g.ch;
    int lun = ppa->g.lun;
    int pl = ppa->g.pl;
    int blk = ppa->g.blk;
    int pg = ppa->g.pg;
    int sec = ppa->g.sec;

    if (ch >= 0 && ch < spp->nchs && lun >= 0 && lun < spp->luns_per_ch && pl >=
        0 && pl < spp->pls_per_lun && blk >= 0 && blk < spp->blks_per_pl && pg
        >= 0 && pg < spp->pgs_per_blk && sec >= 0 && sec < spp->secs_per_pg)
        return true;

    return false;
}

static inline bool valid_lpn(struct ssd *ssd, uint64_t lpn)
{
    return (lpn < ssd->sp.tt_pgs);
}

static inline bool mapped_ppa(struct ppa *ppa)
{
    return !(ppa->ppa == UNMAPPED_PPA);
}

static inline struct ssd_channel *get_ch(struct ssd *ssd, struct ppa *ppa)
{
    return &(ssd->ch[ppa->g.ch]);
}

static inline struct nand_lun *get_lun(struct ssd *ssd, struct ppa *ppa)
{
    struct ssd_channel *ch = get_ch(ssd, ppa);
    return &(ch->lun[ppa->g.lun]);
}

static inline struct nand_plane *get_pl(struct ssd *ssd, struct ppa *ppa)
{
    struct nand_lun *lun = get_lun(ssd, ppa);
    return &(lun->pl[ppa->g.pl]);
}

static inline struct nand_block *get_blk(struct ssd *ssd, struct ppa *ppa)
{
    struct nand_plane *pl = get_pl(ssd, ppa);
    return &(pl->blk[ppa->g.blk]);
}

static inline struct line *get_line(struct ssd *ssd, struct ppa *ppa)
{
    return &(ssd->lm.lines[ppa->g.blk]);
}

static inline struct nand_page *get_pg(struct ssd *ssd, struct ppa *ppa)
{
    struct nand_block *blk = get_blk(ssd, ppa);
    return &(blk->pg[ppa->g.pg]);
}

static uint64_t ssd_advance_status(struct ssd *ssd, struct ppa *ppa, struct
        nand_cmd *ncmd)
{
    int c = ncmd->cmd;
    uint64_t cmd_stime = (ncmd->stime == 0) ? \
        qemu_clock_get_ns(QEMU_CLOCK_REALTIME) : ncmd->stime;
    uint64_t nand_stime;
    struct ssdparams *spp = &ssd->sp;
    struct nand_lun *lun = get_lun(ssd, ppa);
    uint64_t lat = 0;

    switch (c) {
    case NAND_READ:
        /* read: perform NAND cmd first */
        nand_stime = (lun->next_lun_avail_time < cmd_stime) ? cmd_stime : \
                     lun->next_lun_avail_time;
        lun->next_lun_avail_time = nand_stime + spp->pg_rd_lat;
        lat = lun->next_lun_avail_time - cmd_stime;
#if 0
        lun->next_lun_avail_time = nand_stime + spp->pg_rd_lat;

        /* read: then data transfer through channel */
        chnl_stime = (ch->next_ch_avail_time < lun->next_lun_avail_time) ? \
            lun->next_lun_avail_time : ch->next_ch_avail_time;
        ch->next_ch_avail_time = chnl_stime + spp->ch_xfer_lat;

        lat = ch->next_ch_avail_time - cmd_stime;
#endif
        break;

    case NAND_WRITE:
        /* write: transfer data through channel first */
        nand_stime = (lun->next_lun_avail_time < cmd_stime) ? cmd_stime : \
                     lun->next_lun_avail_time;
        if (ncmd->type == USER_IO) {
            lun->next_lun_avail_time = nand_stime + spp->pg_wr_lat;
        } else {
            lun->next_lun_avail_time = nand_stime + spp->pg_wr_lat;
        }
        lat = lun->next_lun_avail_time - cmd_stime;

#if 0
        chnl_stime = (ch->next_ch_avail_time < cmd_stime) ? cmd_stime : \
                     ch->next_ch_avail_time;
        ch->next_ch_avail_time = chnl_stime + spp->ch_xfer_lat;

        /* write: then do NAND program */
        nand_stime = (lun->next_lun_avail_time < ch->next_ch_avail_time) ? \
            ch->next_ch_avail_time : lun->next_lun_avail_time;
        lun->next_lun_avail_time = nand_stime + spp->pg_wr_lat;

        lat = lun->next_lun_avail_time - cmd_stime;
#endif
        break;

    case NAND_ERASE:
        /* erase: only need to advance NAND status */
        nand_stime = (lun->next_lun_avail_time < cmd_stime) ? cmd_stime : \
                     lun->next_lun_avail_time;
        lun->next_lun_avail_time = nand_stime + spp->blk_er_lat;

        lat = lun->next_lun_avail_time - cmd_stime;
        break;

    default:
        ftl_err("Unsupported NAND command: 0x%x\n", c);
    }

    return lat;
}

/* update SSD status about one page from PG_VALID -> PG_VALID */
static void mark_page_invalid(struct ssd *ssd, struct ppa *ppa)
{
    struct line_mgmt *lm = &ssd->lm;
    struct ssdparams *spp = &ssd->sp;
    struct nand_block *blk = NULL;
    struct nand_page *pg = NULL;
    bool was_full_line = false;
    struct line *line;

    /* update corresponding page status */
    pg = get_pg(ssd, ppa);
    ftl_assert(pg->status == PG_VALID);
    pg->status = PG_INVALID;

    /* update corresponding block status */
    blk = get_blk(ssd, ppa);
    ftl_assert(blk->ipc >= 0 && blk->ipc < spp->pgs_per_blk);
    blk->ipc++;
    ftl_assert(blk->vpc > 0 && blk->vpc <= spp->pgs_per_blk);
    blk->vpc--;

    /* update corresponding line status */
    line = get_line(ssd, ppa);
    ftl_assert(line->ipc >= 0 && line->ipc < spp->pgs_per_line);
    if (line->vpc == spp->pgs_per_line) {
        ftl_assert(line->ipc == 0);
        was_full_line = true;
    }
    line->ipc++;
    ftl_assert(line->vpc > 0 && line->vpc <= spp->pgs_per_line);
    /* Adjust the position of the victime line in the pq under over-writes */
    if (line->pos) {
        /* Note that line->vpc will be updated by this call */
        pqueue_change_priority(lm->victim_line_pq, line->vpc - 1, line);
    } else {
        line->vpc--;
    }

    if (was_full_line) {
        /* move line: "full" -> "victim" */
        QTAILQ_REMOVE(&lm->full_line_list, line, entry);
        lm->full_line_cnt--;
        pqueue_insert(lm->victim_line_pq, line);
        lm->victim_line_cnt++;
    }
}

static void mark_page_valid(struct ssd *ssd, struct ppa *ppa)
{
    struct nand_block *blk = NULL;
    struct nand_page *pg = NULL;
    struct line *line;

    /* update page status */
    pg = get_pg(ssd, ppa);
    ftl_assert(pg->status == PG_FREE);
    pg->status = PG_VALID;

    /* update corresponding block status */
    blk = get_blk(ssd, ppa);
    ftl_assert(blk->vpc >= 0 && blk->vpc < ssd->sp.pgs_per_blk);
    blk->vpc++;

    /* update corresponding line status */
    line = get_line(ssd, ppa);
    ftl_assert(line->vpc >= 0 && line->vpc < ssd->sp.pgs_per_line);
    line->vpc++;
}

static void mark_block_free(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    struct nand_block *blk = get_blk(ssd, ppa);
    struct nand_page *pg = NULL;

    for (int i = 0; i < spp->pgs_per_blk; i++) {
        /* reset page status */
        pg = &blk->pg[i];
        ftl_assert(pg->nsecs == spp->secs_per_pg);
        pg->status = PG_FREE;
    }

    /* reset block status */
    ftl_assert(blk->npgs == spp->pgs_per_blk);
    blk->ipc = 0;
    blk->vpc = 0;
    blk->erase_cnt++;
}

static void gc_read_page(struct ssd *ssd, struct ppa *ppa)
{
    /* advance ssd status, we don't care about how long it takes */
    if (ssd->sp.enable_gc_delay) {
        struct nand_cmd gcr;
        gcr.type = GC_IO;
        gcr.cmd = NAND_READ;
        gcr.stime = 0;
        ssd_advance_status(ssd, ppa, &gcr);
    }
}

/* move valid page data (already in DRAM) from victim line to a new page */
static uint64_t gc_write_page(struct ssd *ssd, struct ppa *old_ppa)
{
    struct ppa new_ppa;
    struct nand_lun *new_lun;
    uint64_t lpn = get_rmap_ent(ssd, old_ppa);

    ftl_assert(valid_lpn(ssd, lpn));
    // For garbage collection, we just assign stream 0
    // since our previous guess of lifetime failed
    new_ppa = get_new_page(ssd, 0);
    /* update maptbl */
    set_maptbl_ent(ssd, lpn, &new_ppa);
    /* update rmap */
    set_rmap_ent(ssd, lpn, &new_ppa);

    mark_page_valid(ssd, &new_ppa);

    /* need to advance the write pointer here */
    ssd_advance_write_pointer(ssd, 0);

    if (ssd->sp.enable_gc_delay) {
        struct nand_cmd gcw;
        gcw.type = GC_IO;
        gcw.cmd = NAND_WRITE;
        gcw.stime = 0;
        ssd_advance_status(ssd, &new_ppa, &gcw);
    }

    /* advance per-ch gc_endtime as well */
#if 0
    new_ch = get_ch(ssd, &new_ppa);
    new_ch->gc_endtime = new_ch->next_ch_avail_time;
#endif

    new_lun = get_lun(ssd, &new_ppa);
    new_lun->gc_endtime = new_lun->next_lun_avail_time;

    return 0;
}

static struct line *select_victim_line(struct ssd *ssd, bool force)
{
    struct line_mgmt *lm = &ssd->lm;
    struct line *victim_line = NULL;

    victim_line = pqueue_peek(lm->victim_line_pq);
    if (!victim_line) {
        return NULL;
    }

    if (!force && victim_line->ipc < ssd->sp.pgs_per_line / 8) {
        return NULL;
    }

    pqueue_pop(lm->victim_line_pq);
    victim_line->pos = 0;
    lm->victim_line_cnt--;

    /* victim_line is a danggling node now */
    return victim_line;
}

/* here ppa identifies the block we want to clean */
static void clean_one_block(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    struct nand_page *pg_iter = NULL;
    int cnt = 0;

    for (int pg = 0; pg < spp->pgs_per_blk; pg++) {
        ppa->g.pg = pg;
        pg_iter = get_pg(ssd, ppa);
        /* there shouldn't be any free page in victim blocks */
        ftl_assert(pg_iter->status != PG_FREE);
        if (pg_iter->status == PG_VALID) {
            gc_read_page(ssd, ppa);
            /* delay the maptbl update until "write" happens */
            gc_write_page(ssd, ppa);
            cnt++;
        }
    }

    ssd->pages_from_gc += cnt;

    ftl_assert(get_blk(ssd, ppa)->vpc == cnt);
}

static void mark_line_free(struct ssd *ssd, struct ppa *ppa)
{
    struct line_mgmt *lm = &ssd->lm;
    struct line *line = get_line(ssd, ppa);
    line->ipc = 0;
    line->vpc = 0;
    /* move this line to free line list */
    QTAILQ_INSERT_TAIL(&lm->free_line_list, line, entry);
    lm->free_line_cnt++;
}

static int do_gc(struct ssd *ssd, bool force)
{
    struct line *victim_line = NULL;
    struct ssdparams *spp = &ssd->sp;
    struct nand_lun *lunp;
    struct ppa ppa;
    int ch, lun;

    victim_line = select_victim_line(ssd, force);
    if (!victim_line) {
        return -1;
    }

    ppa.g.blk = victim_line->id;
    ftl_debug("GC-ing line:%d,ipc=%d,victim=%d,full=%d,free=%d\n", ppa.g.blk,
              victim_line->ipc, ssd->lm.victim_line_cnt, ssd->lm.full_line_cnt,
              ssd->lm.free_line_cnt);
    //write_log("GC-ing line:%d,ipc=%d,victim=%d,full=%d,free=%d\n", ppa.g.blk, victim_line->ipc, ssd->lm.victim_line_cnt, ssd->lm.full_line_cnt, ssd->lm.free_line_cnt);

    /* copy back valid data */
    for (ch = 0; ch < spp->nchs; ch++) {
        for (lun = 0; lun < spp->luns_per_ch; lun++) {
            ppa.g.ch = ch;
            ppa.g.lun = lun;
            ppa.g.pl = 0;
            lunp = get_lun(ssd, &ppa);
            clean_one_block(ssd, &ppa);
            mark_block_free(ssd, &ppa);

            if (spp->enable_gc_delay) {
                struct nand_cmd gce;
                gce.type = GC_IO;
                gce.cmd = NAND_ERASE;
                gce.stime = 0;
                ssd_advance_status(ssd, &ppa, &gce);
            }

            lunp->gc_endtime = lunp->next_lun_avail_time;
        }
    }

    /* update line status */
    mark_line_free(ssd, &ppa);

    return 0;
}

static uint64_t ssd_read(struct ssd *ssd, NvmeRequest *req)
{
    struct ssdparams *spp = &ssd->sp;
    uint64_t lba = req->slba;
    int nsecs = req->nlb;
    struct ppa ppa;
    uint64_t start_lpn = lba / spp->secs_per_pg;
    uint64_t end_lpn = (lba + nsecs - 1) / spp->secs_per_pg;
    uint64_t lpn;
    uint64_t sublat, maxlat = 0;

    if (end_lpn >= spp->tt_pgs) {
        ftl_err("start_lpn=%"PRIu64",tt_pgs=%d\n", start_lpn, ssd->sp.tt_pgs);
    }

    /* normal IO read path */
    for (lpn = start_lpn; lpn <= end_lpn; lpn++) {
        ppa = get_maptbl_ent(ssd, lpn);
        if (!mapped_ppa(&ppa) || !valid_ppa(ssd, &ppa)) {
            //printf("%s,lpn(%" PRId64 ") not mapped to valid ppa\n", ssd->ssdname, lpn);
            //printf("Invalid ppa,ch:%d,lun:%d,blk:%d,pl:%d,pg:%d,sec:%d\n",
            //ppa.g.ch, ppa.g.lun, ppa.g.blk, ppa.g.pl, ppa.g.pg, ppa.g.sec);
            continue;
        }

        struct nand_cmd srd;
        srd.type = USER_IO;
        srd.cmd = NAND_READ;
        srd.stime = req->stime;
        sublat = ssd_advance_status(ssd, &ppa, &srd);
        maxlat = (sublat > maxlat) ? sublat : maxlat;
    }

    return maxlat;
}


static uint64_t ssd_write(FemuCtrl *n, struct ssd *ssd, NvmeRequest *req)
{
    NvmeRwCmd *rw = (NvmeRwCmd *) &req->cmd;
    NvmeNamespace *ns = req->ns;
    uint16_t control = le16_to_cpu(rw->control);
    uint32_t dsmgmt = le32_to_cpu(rw->dsmgmt);
    bool stream = control & NVME_RW_DTYPE_STREAMS;
    uint64_t page_death_time = 0;
    uint64_t passed_epoch_since_start = get_passed_epoch_since_start(ssd);
    struct stream_info *si;
    struct stream_info *cmp_si;
    uint16_t dspec = (dsmgmt >> 16) & 0xFFFF; //Stream ID


    // Remember, if stream is true, then dspec = 0 means stream ID 1 (though Linux kernel starts from stream ID = 2),
    // and so on. See nvme_assign_write_stream() line 702 for v5.11.10.
    if (stream){
        // dspec 0 / stream 1 is reserved for data without setting stream ID (WRITE_LIFE_NOT_SET).
        //   or stream 2, for data without a stream ID (WRITE_LIFE_NONE) as defined by Linux Kernel.
        
        // Also, dspec will decrease by 1 for other streams before being passed to the device.
        nvme_update_str_stat(n, ns, dspec);
        ////write_log("NVME_RW_DTYPE_STREAMS Got this for stream: %d\n", dspec);
    }else{
        ////write_log("Stream not set: %d\n", dspec);
        dspec = 0;
    }

    uint64_t lba = req->slba;
    struct ssdparams *spp = &ssd->sp;
    int len = req->nlb;
    uint64_t start_lpn = lba / spp->secs_per_pg;
    uint64_t end_lpn = (lba + len - 1) / spp->secs_per_pg;
    // Add WA tracker
    ssd->pages_from_host += (end_lpn - start_lpn) + 1;
    struct ppa ppa;
    uint64_t lpn;
    uint64_t curlat = 0, maxlat = 0;
    uint64_t chunk;
    uint64_t prediction;
    int stream_choice = 0;
    int i;
    int r;

    //write_log("++++Write start LPN: %"PRIu64", end LPN: %"PRIu64", given stream: %d, now start++++\n", start_lpn, end_lpn, dspec);

    ////write_log("%s, opcode:%#x, start_sec:%#lx, size:%#lx, streamenabled:%d, dspec:%#x\n", __func__, rw->opcode, start_lpn * ssd->sp.secs_per_pg, (end_lpn - start_lpn + 1) * ssd->sp.secsz * ssd->sp.secs_per_pg, stream, dspec);

    if (end_lpn >= spp->tt_pgs) {
        ftl_err("start_lpn=%"PRIu64",tt_pgs=%d\n", start_lpn, ssd->sp.tt_pgs);
    }
 
    while (should_gc_high(ssd)) {
        /* perform GC here until !should_gc(ssd) */
        r = do_gc(ssd, true);
        if (r == -1)
            break;
    }

    // DZ Start
    // This is a write. Update death time and average.
    if (spp->death_time_prediction){
        set_latest_access_time(n, ssd, start_lpn, end_lpn, WRITE_OP);
    }
    // DZ End

    for (lpn = start_lpn; lpn <= end_lpn; lpn++) {
        ppa = get_maptbl_ent(ssd, lpn);
        if (mapped_ppa(&ppa)) {
            /* update old page information first */
            mark_page_invalid(ssd, &ppa);
            set_rmap_ent(ssd, INVALID_LPN, &ppa);
        }
        
        if (dspec > 0){
            stream_choice = dspec;
        }else if (spp->death_time_prediction){
            stream_choice = 0;
            chunk = lpn_to_chunk(lpn, n->pages_per_chunk);
            if (ssd->death_time_list[chunk].last_access_op != INITIAL_OP && ssd->death_time_list[chunk].last_access_op != WRITE_ONCE_OP){
                prediction = ssd->death_time_list[chunk].death_time_avg;
                page_death_time = passed_epoch_since_start + prediction;
                // Stream should have exponential accepted range
                // e.g. stream 1: DT b/t 0~1
                //      stream 2: DT b/t 1~3
                //      stream 3: DT b/t 3~7
                for (i = 0; i < spp->msl; i++){
                    prediction >>= 1;
                    stream_choice += 1;
                    if (prediction == 0){
                        break;
                    }
                }

                // If lifetime is too long, and we cannot handle them,
                // then give this to default stream
                if (prediction > 0){
                    stream_choice = 0;
                }
                si = &ssd->stream_info[stream_choice];
                if (spp->enable_stream_redirect){
                    prediction = ssd->death_time_list[chunk].death_time_avg;
                    if (stream_choice > 0 && si->receiver == false){
                        // Check if there is another write pointer we can redirect this page to
                        // Also, if this stream has received something from some other stream, then we don't redirect anything from this stream to other streams
                        // AKA receiver should not send redirects
                        for (i = 1; i <= spp->msl; i++){
                            cmp_si = &ssd->stream_info[i];
                            // Sender should not receive
                            if (cmp_si->sender){
                                continue;
                            }
                            // Check if L > (P - 1) * V_i
                            if ((pow(2, i) - 1) * spp->access_interval_precision > (spp->pages_per_superblock - 1) * cmp_si->avg_incoming_interval){
                                // The target must have L > (P - 1) * V_i, goes here
                                if (page_death_time >= si->earliest_death_time && page_death_time <= si->latest_death_time){
                                    // Redirect
                                    write_log("Page lifetime prediction: %"PRIu64", deathtime prediction: %"PRIu64",\nredirected from stream %d (T_e = %"PRIu64", T_l = %"PRIu64", T_o = %"PRIu64") to %d (T_e = %"PRIu64", T_l = %"PRIu64", T_o = %"PRIu64").\n\n", prediction, page_death_time, stream_choice, si->earliest_death_time, si->latest_death_time, si->block_start_time, i, cmp_si->earliest_death_time, cmp_si->latest_death_time, cmp_si->block_start_time);
                                    stream_choice = i;
                                    si->sender = true;
                                    cmp_si->receiver = true;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                //write_log("Addr: %"PRIu64", Chunk: %"PRIu64", Avg Life: %d, Assigned to stream: %d.\n", lpn, chunk, ssd->death_time_list[chunk].death_time_avg, stream_choice);
            }else if(ssd->death_time_list[chunk].last_access_op == WRITE_ONCE_OP){
                stream_choice = 0;
                //write_log("Addr: %"PRIu64", Chunk: %"PRIu64", First time written, no DT info, assigned to stream 0.\n", lpn, chunk);
            }
        }
        si = &ssd->stream_info[stream_choice];
        // Change the write pointer earliest/latest death time
        if (page_death_time < si->earliest_death_time || si->page_counter == 0){
            si->earliest_death_time = page_death_time;
        }
        if (page_death_time > si->latest_death_time){
            si->latest_death_time = page_death_time;
        }

        /* new write */
        ppa = get_new_page(ssd, stream_choice);
        /* update maptbl */
        set_maptbl_ent(ssd, lpn, &ppa);
        /* update rmap */
        set_rmap_ent(ssd, lpn, &ppa);

        mark_page_valid(ssd, &ppa);

        /* need to advance the write pointer here */
        ssd_advance_write_pointer(ssd, stream_choice);

        struct nand_cmd swr;
        swr.type = USER_IO;
        swr.cmd = NAND_WRITE;
        swr.stime = req->stime;
        /* get latency statistics */
        curlat = ssd_advance_status(ssd, &ppa, &swr);
        maxlat = (curlat > maxlat) ? curlat : maxlat;
    }
    //write_log("----Write start LPN: %"PRIu64", end LPN: %"PRIu64", given stream: %d, now end----\n", start_lpn, end_lpn, dspec);
    return maxlat;
}

// DZ Start
static void ssd_dsm(FemuCtrl *n, struct ssd *ssd, NvmeRequest *req){
    // Mostly adapted from nvme_dsm() in nvme-io.c.
    uint64_t lba = req->slba;
    struct ssdparams *spp = &ssd->sp;
    int len = req->nlb;
    struct ppa ppa;
    uint64_t lpn;
    uint64_t start_lpn;
    uint64_t end_lpn;

    if (req->cmd.cdw11 & NVME_DSMGMT_AD){
        uint16_t nr = (req->cmd.cdw10 & 0xff) + 1;
        NvmeDsmRange range[nr];
        
        // FEMU will handle the real I/O request first
        // and also finished all sanity check on DSM range.
        // See nvme_dsm() in nvme-io.c.
        // However, we still need to get range information using this function.
        dma_write_prp(n, (uint8_t *)range, sizeof(range), req->cmd.dptr.prp1, req->cmd.dptr.prp2);
        // We can skip sanity check here.
        for (int i = 0; i < nr; i++) {
            lba = le64_to_cpu(range[i].slba);
            len = le32_to_cpu(range[i].nlb);

            start_lpn = lba / spp->secs_per_pg;
            end_lpn = (lba + len - 1) / spp->secs_per_pg;
            if (end_lpn >= spp->tt_pgs) {
                ftl_err("start_lpn=%"PRIu64",tt_pgs=%d\n", start_lpn, ssd->sp.tt_pgs);
            }
            // Mark these pages as invalid
            for (lpn = start_lpn; lpn <= end_lpn; lpn++) {
                ppa = get_maptbl_ent(ssd, lpn);
                if (mapped_ppa(&ppa)) {
                    // This physical page is invalid now
                    mark_page_invalid(ssd, &ppa);
                    // This LPN is also invalid now
                    ssd->maptbl[lpn].ppa = UNMAPPED_PPA;
                    set_rmap_ent(ssd, INVALID_LPN, &ppa);
                }
            }
            if (spp->death_time_prediction){
                set_latest_access_time(n, ssd, start_lpn, end_lpn, DISCARD_OP);
            }
        }
    }
}
// DZ End

static void *ftl_thread(void *arg)
{
    #ifdef FEMU_DEBUG_FTL
    femu_log_file = fopen("/tmp/femu.log","a");
    #endif
    FemuCtrl *n = (FemuCtrl *)arg;
    struct ssd *ssd = n->ssd;
    NvmeRequest *req = NULL;
    uint64_t lat = 0;
    int rc;
    int i;

    while (!*(ssd->dataplane_started_ptr)) {
        usleep(100000);
    }

    /* FIXME: not safe, to handle ->to_ftl and ->to_poller gracefully */
    ssd->to_ftl = n->to_ftl;
    ssd->to_poller = n->to_poller;

    while (1) {
        for (i = 1; i <= n->num_poller; i++) {
            if (!ssd->to_ftl[i] || !femu_ring_count(ssd->to_ftl[i]))
                continue;

            rc = femu_ring_dequeue(ssd->to_ftl[i], (void *)&req, 1);
            if (rc != 1) {
                printf("FEMU: FTL to_ftl dequeue failed\n");
            }

            ftl_assert(req);
            switch (req->cmd.opcode) {
            case NVME_CMD_WRITE:
                lat = ssd_write(n, ssd, req);
                break;
            case NVME_CMD_READ:
                lat = ssd_read(ssd, req);
                break;
            case NVME_CMD_DSM:
                lat = 0;
                // DZ Start
                // Handle discard request here
                ssd_dsm(n, ssd, req);
                // DZ End
                break;
            default:
                //ftl_err("FTL received unkown request type, ERROR\n");
                ;
            }

            req->reqlat = lat;
            req->expire_time += lat;

            rc = femu_ring_enqueue(ssd->to_poller[i], (void *)&req, 1);
            if (rc != 1) {
                ftl_err("FTL to_poller enqueue failed\n");
            }

            /* clean one line if needed (in the background) */
            if (should_gc(ssd)) {
                do_gc(ssd, false);
            }
        }
    }
    #ifdef FEMU_DEBUG_FTL
    fclose(femu_log_file);
    #endif
    return NULL;
}

