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

#define max(a, b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a, b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

static void *ftl_thread(void *arg);
static inline struct nand_block *get_blk(struct ssd *ssd, struct ppa *ppa);
// static void pqueue_change_line_lun(struct line * line, int curchannel, int curlun)
// {
//     line->curchannel = curchannel;
//     line->curlun = curlun;
// }

static int get_channels_needed(struct ssd *ssd, int stream){
    struct ssdparams *spp = &ssd->sp;
    return spp->default_channels_per_line;
}

static int get_luns_needed(struct ssd *ssd, int stream){
    struct ssdparams *spp = &ssd->sp;
    if (spp->enable_hetero_sbsize){
        if (stream == spp->gc_stream_id){
            return spp->luns_per_ch;
        }
        if (stream % 2 == 0){
            return 1;
        }else{
            return spp->default_luns_per_channel;
        }
    }
    return spp->default_luns_per_channel;
}

static bool should_gc_channel_lun(struct ssd *ssd, int channel, int lun, bool high){
    struct ssdparams *spp = &ssd->sp;
    if (ssd->channel_mgmt.channel[channel].lun[lun].free_blocks_cnt < (high ? spp->gc_thres_blocks_high : spp->gc_thres_blocks)){
        return true;
    }
    return false;
}

static bool should_gc(struct ssd *ssd)
{
    int i, j;
    struct ssdparams *spp = &ssd->sp;
    for (i = 0; i < spp->nchs; i++){
        for (j = 0; j < spp->luns_per_ch; j++){
            // write_log("channel: %d, lun: %d, free_blocks: %d, threshold: %d\n", i, j, ssd->channel_mgmt.channel[i].lun[j].free_blocks_cnt, spp->gc_thres_blocks);
            if (should_gc_channel_lun(ssd, i, j, false)){
                // write_log("should gc channel %d, lun %d\n", i, j);
                return true;
            }
        }
    }
    return false;
}

static bool should_gc_high(struct ssd *ssd)
{
    int i, j;
    struct ssdparams *spp = &ssd->sp;
    for (i = 0; i < spp->nchs; i++){
        for (j = 0; j < spp->luns_per_ch; j++){
            if (should_gc_channel_lun(ssd, i, j, true)){
                return true;
            }
        }
    }
    return false;
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

static struct seq_write_info *init_seq_write_info(const uint size){
    struct seq_write_info *s = g_malloc0(sizeof(struct seq_write_info));
    s->size = size;
    s->cur = 0;
    s->list = g_malloc0_n(size, sizeof(uint64_t));
    s->inited = 0;
    return s;
}

static bool check_seq_write_info(struct seq_write_info *s, uint64_t addr){
    for (int i = 0; i < s->size; i++){
        if (unlikely(i == s->inited)){
            s->list[i] = addr;
            s->inited += 1;
            s->cur += 1;
            return false;
        }else{
            if (s->list[i] == addr){
                s->list[i] = addr;
                return true;
            }
        }
    }
    s->cur %= s->size;
    s->list[s->cur] = addr;
    s->cur += 1;
    return false;
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
    uint64_t *chunks = g_malloc0_n(((end_lpn - start_lpn) + 1), sizeof(uint64_t));
    uint64_t num_chunks = get_chunk_list(n, ssd, start_lpn, end_lpn, chunks);
    bool is_recent;
    for (i = 0; i < num_chunks; i++){
        chunk = chunks[i];
        is_recent = check_seq_write_info(ssd->seq_info, chunk);
        if (is_recent && op == WRITE_OP){
            continue;
        }
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

// static void victim_line_print(FILE *out, void *a){
//     struct line * line = (struct line *)a;
//     write_log("pos: %"PRId64", vpc: %d\n", line->pos[line->curchannel][line->curlun], line->vpc);
// }

// static size_t victim_line_get_pos(void *a)
// {
//     struct line * line = (struct line *)a;
//     return line->pos[line->curchannel][line->curlun];
// }

// static void victim_line_set_pos(void *a, size_t pos)
// {
//     struct line * line = (struct line *)a;
//     line->pos[line->curchannel][line->curlun] = pos;
// }

static struct line* create_line(struct ssd *ssd, int stream){
    struct ssdparams *spp = &ssd->sp;
    struct line_mgmt *lm = &ssd->lm;
    int i, j, k;
    int channels = get_channels_needed(ssd, stream);
    int luns = get_luns_needed(ssd, stream);
    struct line *line = (struct line *)g_malloc0(sizeof(struct line));
    struct channel_mgmt *channel_mgmt = &ssd->channel_mgmt;
    struct lun_mgmt *lun_mgmt;
    struct block_mgmt *block_mgmt;
    struct block_num *block_num;
    struct ppa ppa;
    line->id = channel_mgmt->next_line_id;
    write_log("[7, %d, %d, %d]\n", channels, luns, line->id);
    line->ipc = 0;
    line->vpc = 0;
    line->valid = true;
    line->use = USE_INUSE;
    line->stream = stream;
    line->total_channels = channels;
    line->pgs_per_line = 0;
    line->channel_list = g_malloc0_n(line->total_channels, sizeof(int));
    line->total_luns = g_malloc0_n(line->total_channels, sizeof(int));
    line->lun_list = g_malloc0_n(line->total_channels, sizeof(int*));
    line->block_list = g_malloc0_n(line->total_channels, sizeof(int*));
    line->full_entry = g_malloc0_n(spp->nchs, sizeof(QTAILQ_ENTRY(line) *));
    line->victim_entry = g_malloc0_n(spp->nchs, sizeof(QTAILQ_ENTRY(line) *));
    line->inserted_to_full_queue = g_malloc0_n(spp->nchs, sizeof(bool*));
    line->inserted_to_victim_queue = g_malloc0_n(spp->nchs, sizeof(bool*));
    for (i = 0; i < spp->nchs; i++){
        line->full_entry[i] = g_malloc0_n(spp->luns_per_ch, sizeof(QTAILQ_ENTRY(line)));
        line->victim_entry[i] = g_malloc0_n(spp->luns_per_ch, sizeof(QTAILQ_ENTRY(line)));
        line->inserted_to_full_queue[i] = g_malloc0_n(spp->luns_per_ch, sizeof(bool));
        line->inserted_to_victim_queue[i] = g_malloc0_n(spp->luns_per_ch, sizeof(bool));
        for (j = 0; j < spp->luns_per_ch; j++){
            line->inserted_to_full_queue[i][j] = false;
            line->inserted_to_victim_queue[i][j] = false;
        }
    }
    ppa.g.pl = 0;
    for (j = 0; j < line->total_channels; j++){
        // write_log("j = %d\n", j);
        
        lun_mgmt = &channel_mgmt->channel[channel_mgmt->next_avail_channel];
        // write_log("next_avail_channel = %d\n", channel_mgmt->next_avail_channel);
        
        // Assign flexible channels
        line->channel_list[j] = channel_mgmt->next_avail_channel;
        // Assign flexible LUNs
        line->total_luns[j] = luns;
        
        line->lun_list[j] = g_malloc0_n(line->total_luns[j], sizeof(int));
        line->block_list[j] = g_malloc0_n(line->total_luns[j], sizeof(int));
        for (k = 0; k < line->total_luns[j]; k++){
            block_mgmt = &lun_mgmt->lun[lun_mgmt->next_avail_lun];
            if (block_mgmt->free_blocks_cnt == 0){
                continue;
            }
            block_num = QTAILQ_FIRST(&block_mgmt->free_block_list);
            QTAILQ_REMOVE(&block_mgmt->free_block_list, block_num, entry);
            block_mgmt->free_blocks_cnt--;
            line->lun_list[j][k] = lun_mgmt->next_avail_lun;
            line->block_list[j][k] = block_num->block_num;
            g_free(block_num);
            ppa.g.ch = line->channel_list[j];
            ppa.g.lun = line->lun_list[j][k];
            ppa.g.blk = line->block_list[j][k];
            get_blk(ssd, &ppa)->status = USE_INUSE;
            lm->channel_lines[line->block_list[j][k]][line->channel_list[j]][line->lun_list[j][k]] = line;
            line->pgs_per_line += spp->pgs_per_blk;
            lun_mgmt->next_avail_lun = (lun_mgmt->next_avail_lun + 1) % spp->luns_per_ch;
        }
        channel_mgmt->next_avail_channel = (channel_mgmt->next_avail_channel + 1) % spp->nchs;
    }
    channel_mgmt->next_line_id++;
    return line;
}

static void ssd_init_lines(struct ssd *ssd)
{
    write_log("ssd_init_lines start\n");
    
    struct ssdparams *spp = &ssd->sp;
    struct line_mgmt *lm = &ssd->lm;
    struct channel_mgmt *channel_mgmt = &ssd->channel_mgmt;
    struct lun_mgmt *lun_mgmt;
    struct block_mgmt *block_mgmt;
    struct block_num *block_num;
    // spp->max_allow_gc_lines = 0;
    channel_mgmt->next_avail_channel = 0;
    channel_mgmt->next_line_id = 0;
    channel_mgmt->channel = g_malloc0_n(spp->nchs, sizeof(struct lun_mgmt));
    // struct line *line;
    int i, j, k;
    lm->line_resource = g_malloc0_n(spp->nchs, sizeof(struct line_resource_mgmt*));
    for (i = 0; i < spp->nchs; i++){
        lm->line_resource[i] = g_malloc0_n(spp->luns_per_ch, sizeof(struct line_resource_mgmt));
        lun_mgmt = &channel_mgmt->channel[i];
        lun_mgmt->next_avail_lun = 0;
        lun_mgmt->lun = g_malloc0_n(spp->luns_per_ch, sizeof(struct block_mgmt));
        for (j = 0; j < spp->luns_per_ch; j++){
            QTAILQ_INIT(&lm->line_resource[i][j].free_line_list);
            QTAILQ_INIT(&lm->line_resource[i][j].full_line_list);
            QTAILQ_INIT(&lm->line_resource[i][j].victim_line_list);
            lm->line_resource[i][j].free_line_cnt = 0;
            lm->line_resource[i][j].victim_line_cnt = 0;
            lm->line_resource[i][j].full_line_cnt = 0;
            block_mgmt = &lun_mgmt->lun[j];
            block_mgmt->free_blocks_cnt = spp->blks_per_lun;
            QTAILQ_INIT(&block_mgmt->free_block_list);
            for (k = 0; k < spp->blks_per_pl; k++){
                block_num = g_malloc0(sizeof(struct block_num));
                block_num->block_num = k;
                QTAILQ_INSERT_TAIL(&block_mgmt->free_block_list, block_num, entry);
            }
        }
    }
    lm->tt_lines = spp->tt_lines;
    ftl_assert(lm->tt_lines == spp->tt_lines);
    //lm->lines = g_malloc0(sizeof(struct line) * lm->tt_lines);
    lm->channel_lines = g_malloc0_n(spp->blks_per_lun, sizeof(struct line***));
    for (i = 0; i < spp->blks_per_lun; i++){
        lm->channel_lines[i] = g_malloc0_n(spp->nchs, sizeof(struct line**));
        for (j = 0; j < spp->nchs; j++){
            for (k = 0; k < spp->luns_per_ch; k++){
                lm->channel_lines[i][j] = g_malloc0_n(spp->luns_per_ch, sizeof(struct line*));
            }
        }
    }
    write_log("ssd_init_lines end\n");
    
}

static void ssd_init_write_pointer(struct ssd *ssd, uint8_t streams)
{
    write_log("ssd_init_write_pointer start\n");
    
    struct ssdparams *spp = &ssd->sp;
    spp->real_num_streams = streams + 2;
    if (spp->enable_hetero_sbsize){
        spp->real_num_streams = streams + 2 + 1; // 1 for gc, 2 for the default stream (8/1, 8/8 or maybe 8/4)
    }
    
    spp->gc_stream_id = spp->real_num_streams - 1;
    // n streams -> n+2 write pointers since we need stream 0 as default stream and stream n+1 as GC stream.
    // * 2 for shadow streams 
    ssd->wp = g_malloc0_n(spp->real_num_streams, sizeof(struct write_pointer));
    ssd->stream_info = g_malloc0_n(spp->real_num_streams, sizeof(struct stream_info));
    ssd->seq_info = init_seq_write_info(streams * 4);
    struct write_pointer *wpp = NULL;
    struct stream_info *si = NULL;
    // struct line_mgmt *lm = &ssd->lm;
    struct line *curline = NULL;
    for (int i = 0; i < spp->real_num_streams; i++){
        //curline = QTAILQ_FIRST(&lm->free_line_list);
        curline = create_line(ssd, i);
        // QTAILQ_REMOVE(&lm->free_line_list, curline, entry);
        // lm->free_line_cnt--;
        wpp = &ssd->wp[i];
        si = &ssd->stream_info[i];
        /* wpp->curline is always our next-to-write super-block */
        wpp->curline = curline;
        wpp->ch = 0;
        wpp->lun = 0;
        wpp->pg = 0;
        wpp->blk = 0;
        wpp->pl = 0;
        si->earliest_death_time = -1; // Overflow on purpose. Max uint64_t.
        si->latest_death_time = 0;
        si->avg_incoming_interval = 0;
        si->avg_temp_incoming_interval = 0;
        si->block_open_time = get_uptime(ssd);
        si->stream_counter_start_time = get_uptime(ssd);
        si->full_before = false;
        si->sender = false;
        si->receiver = false;
        si->next_avail_time = qemu_clock_get_ns(QEMU_CLOCK_REALTIME);
    }
    write_log("ssd_init_write_pointer end\n");
    
}

static inline void check_addr(int a, int max)
{
    ftl_assert(a >= 0 && a < max);
}

static struct line *get_next_free_line(struct ssd *ssd, int stream)
{
    // struct line_mgmt *lm = &ssd->lm;
    struct line *curline = NULL;

    curline = create_line(ssd, stream);
    if (!curline) {
        ftl_err("No free lines left in [%s] !!!!\n", ssd->ssdname);
        return NULL;
    }

    return curline;
}

static void ssd_advance_write_pointer(struct ssd *ssd, uint8_t stream, uint64_t lpn)
{
    // write_log("debug 14.1, stream = %d\n", stream);
    
    struct ssdparams *spp = &ssd->sp;
    struct write_pointer *wpp = &ssd->wp[stream];
    struct stream_info *si = &ssd->stream_info[stream];
    uint64_t uptime = get_uptime(ssd);
    struct line_mgmt *lm = &ssd->lm;
    int real_channel_no;
    int real_lun_no;
    int i, j;
    struct ppa ppa;
    ppa.g.pl = 0;
    // write_log("debug 14.2\n");
    
    //check_addr(wpp->ch, wpp->curline->start_channel + wpp->curline->total_channels);
    if (!(wpp->curline->channel_list[wpp->ch] >= 0 && wpp->curline->channel_list[wpp->ch] < spp->nchs)){
        ftl_debug("ch bug: wpp->ch: %d, wpp->curline->channel_list[wpp->ch]: %d, spp->nchs: %d\n", wpp->ch, wpp->curline->channel_list[wpp->ch], spp->nchs);
    }
    check_addr(wpp->curline->channel_list[wpp->ch], spp->nchs);
    wpp->ch++;
    if (wpp->ch == wpp->curline->total_channels) {
        wpp->ch = 0;
        // write_log("debug 14.3\n");
        
        if (!(wpp->lun >= 0 && wpp->lun < wpp->curline->total_luns[wpp->ch])){
            ftl_debug("lun bug: wpp->lun: %d, luns of this channel: %d\n", wpp->lun, wpp->curline->total_luns[wpp->ch]);
        }
        check_addr(wpp->lun, wpp->curline->total_luns[wpp->ch]);
        wpp->lun++;
        /* in this case, we should go to next lun */
        if (wpp->lun == wpp->curline->total_luns[wpp->ch]) {
            wpp->lun = 0;
            // write_log("debug 14.4\n");
            
            /* go to next page in the block */
            if (!(wpp->pg >= 0 && wpp->pg < spp->pgs_per_blk)){
                ftl_debug("pg bug: wpp->pg: %d, spp->pgs_per_blk: %d\n", wpp->pg, spp->pgs_per_blk);
            }
            check_addr(wpp->pg, spp->pgs_per_blk);
            wpp->pg++;
            if (wpp->pg == spp->pgs_per_blk) {
                wpp->pg = 0;
                /* move current line to {victim,full} line list */
                // write_log("debug 14.5\n");
                
                for (i = 0; i < wpp->curline->total_channels; i++){
                    real_channel_no = wpp->curline->channel_list[i];
                    for (j = 0; j < wpp->curline->total_luns[i]; j++){
                        real_lun_no = wpp->curline->lun_list[i][j];
                        ppa.g.ch = real_channel_no;
                        ppa.g.lun = real_lun_no;
                        ppa.g.blk = wpp->curline->block_list[i][j];
                        get_blk(ssd, &ppa)->status = USE_FULL;
                        // write_log("debug 14.6 i = %d/%d, j = %d/%d, channel = %d, lun = %d\n", i, wpp->curline->total_channels, j, wpp->curline->total_luns[i], real_channel_no, real_lun_no);
                        
                        if (wpp->curline->vpc == wpp->curline->pgs_per_line) {
                            /* all pgs are still valid, move to full line list */
                            // write_log("debug 14.7 put in full, wpp->curline->id = %d\n", wpp->curline->id);
                            
                            ftl_assert(wpp->curline->ipc == 0);
                            // victim_line_assert(ssd, real_channel_no, real_lun_no, wpp->curline);
                            if (wpp->curline->inserted_to_full_queue[real_channel_no][real_lun_no]){

                            }else{
                                QTAILQ_INSERT_TAIL(&lm->line_resource[real_channel_no][real_lun_no].full_line_list, wpp->curline, full_entry[real_channel_no][real_lun_no]);
                                wpp->curline->inserted_to_full_queue[real_channel_no][real_lun_no] = true;
                                lm->line_resource[real_channel_no][real_lun_no].full_line_cnt++;
                            }
                            // victim_line_assert(ssd, real_channel_no, real_lun_no, wpp->curline);
                            // write_log("debug 14.8\n");
                            
                            
                        } else {
                            ftl_assert(wpp->curline->vpc >= 0 && wpp->curline->vpc < wpp->curline->pgs_per_line);
                            /* there must be some invalid pages in this line */
                            ftl_assert(wpp->curline->ipc > 0);
                            if (wpp->curline->inserted_to_victim_queue[real_channel_no][real_lun_no]){

                            }else{
                                QTAILQ_INSERT_TAIL(&lm->line_resource[real_channel_no][real_lun_no].victim_line_list, wpp->curline, victim_entry[real_channel_no][real_lun_no]);
                                wpp->curline->inserted_to_victim_queue[real_channel_no][real_lun_no] = true;
                                lm->line_resource[real_channel_no][real_lun_no].victim_line_cnt++;
                            }
                        }
                    }
                }
                /* current line is used up, pick another empty line */

                // DZ Start
                wpp->curline->use = USE_FULL;
                wpp->curline->close_time = uptime;
                wpp->curline->stream = stream;
                wpp->curline->earliest_dt = si->earliest_death_time;
                wpp->curline->latest_dt = si->latest_death_time;
                write_log("Stream block %d closed, earliest DT: %"PRIu64", latest DT: %"PRIu64", Now: %"PRIu64"\n", stream, si->earliest_death_time, si->latest_death_time, uptime);
                //wpp->curline->expected_h = ;
                si->sender = false;
                si->receiver = false;
                si->block_open_time = uptime;
                si->earliest_death_time = -1; // Max uint64_t
                si->latest_death_time = 0;
                si->next_avail_time = qemu_clock_get_ns(QEMU_CLOCK_REALTIME);
                wpp->curline = NULL;
                // write_log("debug 14.13\n");
                
                wpp->curline = get_next_free_line(ssd, stream);
                /* ftl_debug("Stream %d now has line %d,victim=%d,full=%d,free=%d,lpn=%"PRIu64",hostpgs=%"PRIu64",gcpgs=%"PRIu64"\n", stream, wpp->curline->id,
                    ssd->lm.victim_line_cnt, ssd->lm.full_line_cnt,
                    ssd->lm.free_line_cnt, lpn, ssd->pages_from_host, ssd->pages_from_gc);*/
                //write_log("[4, %d, %d, %d, %d, %d, %"PRIu64", %"PRIu64"]\n", stream, wpp->curline->id, ssd->lm.victim_line_cnt, ssd->lm.full_line_cnt, ssd->lm.free_line_cnt, ssd->pages_from_host, ssd->pages_from_gc);
                if (!wpp->curline) {
                    /* TODO */
                    ftl_debug("Failed to get a line\n");
                    fflush(NULL);
                    abort();
                }
                // write_log("debug 14.15\n");
                
                wpp->ch = 0;
                wpp->blk = 0;
                //femu_log("New superblock: %d caused by stream: %d\n", wpp->blk, stream);
                if (!(wpp->blk >= 0 && wpp->blk < spp->blks_per_pl)){
                    ftl_debug("blk bug: wpp->blk: %d, spp->blks_per_pl: %d\n", wpp->blk, spp->blks_per_pl);
                }
                // write_log("debug 14.16\n");
                
                check_addr(wpp->blk, spp->blks_per_pl);
                // write_log("debug 14.17\n");
                
                /* make sure we are starting from page 0 in the super block */
                ftl_assert(wpp->pg == 0);
                ftl_assert(wpp->lun == 0);
                ftl_assert(wpp->ch == 0);
                /* TODO: assume # of pl_per_lun is 1, fix later */
                ftl_assert(wpp->pl == 0);
            }
        }
    }
    // write_log("debug 14.18\n");
    
}

static struct ppa get_new_page(struct ssd *ssd, uint8_t stream)
{
    struct write_pointer *wpp = &ssd->wp[stream];
    struct ppa ppa;
    ppa.ppa = 0;
    ppa.g.ch = wpp->curline->channel_list[wpp->ch];
    ppa.g.lun = wpp->curline->lun_list[wpp->ch][wpp->lun];
    ppa.g.pg = wpp->pg;
    ppa.g.blk = wpp->curline->block_list[wpp->ch][wpp->lun];
    ppa.g.pl = wpp->pl;
    // write_log("wpp->ch: %d, ", wpp->ch);
    // write_log("Stream: %d, channel: %d, lun: %d, pg: %d, blk: %d, Plane: %d\n", stream, ppa.g.ch, ppa.g.lun, ppa.g.pg, ppa.g.blk, ppa.g.pl);
    
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
    // spp->secsz = 512;
    // spp->secs_per_pg = 8;
    // spp->pgs_per_blk = 1024;
    // spp->blks_per_pl = 512; /* 128GB */

    // Old SSD
    // Page size: 4K
    // 64 Pages per block
    // 1024 Blocks
    // Total: 16GB
    // spp->secsz = 512;
    // spp->secs_per_pg = 8;
    // spp->pgs_per_blk = 64;
    // spp->blks_per_pl = 1024;

    // Newer SSD
    // Page size: 8K
    // 128 Pages per block
    // 1024 Blocks
    // Total: 64GB
    // spp->secsz = 512;
    // spp->secs_per_pg = 16;
    // spp->pgs_per_blk = 128;
    // spp->blks_per_pl = 1024;

    // Latest SSD
    // Page size: 16K
    // 256 Pages per block
    // 1024 Blocks
    // Total: 256GB
    spp->secsz = 512;
    spp->secs_per_pg = 32;
    spp->pgs_per_blk = 256;
    spp->blks_per_pl = 1024;

    //spp->blks_per_pl = 320; /* 20GB */
    spp->pls_per_lun = 1;
    spp->luns_per_ch = 8;
    spp->nchs = 8;
    spp->gc_start_channel = 0;
    spp->gc_start_lun = 0;
    //spp->min_channels_per_line = spp->nchs / spp->channel_regions;
    //spp->min_channels_per_line = spp->default_channels_per_line;

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
    ftl_assert(spp->nchs % spp->default_channels_per_line == 0);
    ftl_assert(spp->luns_per_ch % spp->default_luns_per_channel == 0);
    ftl_assert(spp->blks_per_pl % spp->init_blk_per_plane == 0);
    //spp->channel_regions = (1 << spp->channel_split_exp);
    spp->channel_regions = spp->nchs / spp->default_channels_per_line;
    spp->lun_regions = spp->luns_per_ch / spp->default_luns_per_channel;
    spp->blk_regions = spp->blks_per_pl / spp->init_blk_per_plane;
    /* line is special, put it at the end */
    spp->blks_per_line = spp->tt_luns; /* TODO: to fix under multiplanes */
    spp->pgs_per_line = spp->blks_per_line * spp->pgs_per_blk;
    spp->secs_per_line = spp->pgs_per_line * spp->secs_per_pg;
    spp->tt_lines = spp->blks_per_lun * spp->channel_regions * spp->lun_regions / spp->init_blk_per_plane; /* TODO: to fix under multiplanes */

    spp->gc_thres_pcent = 0.75;
    spp->gc_thres_blocks = (int)((1 - spp->gc_thres_pcent) * spp->blks_per_pl);
    spp->gc_thres_pcent_high = 0.95;
    spp->gc_thres_blocks_high = (int)((1 - spp->gc_thres_pcent_high) * spp->blks_per_pl);
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
    blk->status = USE_FREE;
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

    ssd->death_time_list = g_malloc0_n(spp->tt_chunks, sizeof(struct death_time_track));
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
    #ifdef FEMU_DEBUG_FTL
    char str[80];
    sprintf(str, "/mnt/testpartition/femu%d.log", n->virt_id);
    femu_log_file = fopen(str, "w+");
    #endif
    struct ssd *ssd = n->ssd;
    ssd->pages_from_host = 0;
    ssd->pages_from_gc = 0;

    struct ssdparams *spp = &ssd->sp;

    ftl_assert(ssd);
    spp->channel_split_exp = n->channel_split_exp;
    spp->default_channels_per_line = n->default_channels_per_line;
    spp->default_luns_per_channel = n->default_luns_per_channel;
    spp->init_blk_per_plane = n->init_blk_per_plane;
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
    if (!spp->enable_stream){
        spp->msl = 0;
    }
    spp->enable_stream_redirect = n->enable_stream_redirect;
    spp->access_interval_precision = n->access_interval_precision;
    spp->enable_hetero_sbsize = n->enable_hetero_sbsize;
    spp->enable_partial_gc = n->enable_partial_gc;
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
    return (ssd->lm.channel_lines[ppa->g.blk][ppa->g.ch][ppa->g.lun]);
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
    // write_log("debug 8.3.1\n");
    
    struct line_mgmt *lm = &ssd->lm;
    struct ssdparams *spp = &ssd->sp;
    struct nand_block *blk = NULL;
    struct nand_page *pg = NULL;
    int channel, lun;
    bool was_full_line = false;
    struct line *line;
    int i, j;
    // write_log("debug 8.3.2\n");
    
    /* update corresponding page status */
    pg = get_pg(ssd, ppa);
    ftl_assert(pg->status == PG_VALID);
    pg->status = PG_INVALID;
    // write_log("debug 8.3.4\n");
    
    /* update corresponding block status */
    blk = get_blk(ssd, ppa);
    ftl_assert(blk->ipc >= 0 && blk->ipc < spp->pgs_per_blk);
    blk->ipc++;
    ftl_assert(blk->vpc > 0 && blk->vpc <= spp->pgs_per_blk);
    blk->vpc--;
    // write_log("debug 8.3.4\n");
    
    /* update corresponding line status */
    line = get_line(ssd, ppa);
    ftl_assert(line->ipc >= 0 && line->ipc < line->pgs_per_line);
    if (line->vpc == line->pgs_per_line) {
        ftl_assert(line->ipc == 0);
        was_full_line = true;
    }
    // write_log("debug 8.3.5\n");
    
    line->ipc++;
    ftl_assert(line->vpc > 0 && line->vpc <= line->pgs_per_line);
    /* Adjust the position of the victime line in the pq under over-writes */
    line->vpc--;
    // write_log("debug 8.3.6\n");
    
    for (i = 0; i < line->total_channels; i++){
        channel = line->channel_list[i];
        for (j = 0; j < line->total_luns[i]; j++){
            lun = line->lun_list[i][j];
            // write_log("debug 8.3.7, i = %d/%d, j = %d/%d, channel = %d, lun = %d\n", i, line->total_channels, j, line->total_luns[i], channel, lun);
            
            
            if (was_full_line) {
                /* move line: "full" -> "victim" */
                // write_log("debug 8.3.10\n");
                
                if (line->inserted_to_full_queue[channel][lun]){
                    QTAILQ_REMOVE(&lm->line_resource[channel][lun].full_line_list, line, full_entry[channel][lun]);
                    line->inserted_to_full_queue[channel][lun] = false;
                    lm->line_resource[channel][lun].full_line_cnt--;
                }else{

                }
                // write_log("debug 8.3.11\n");
                
                if (line->inserted_to_victim_queue[channel][lun]){

                }else{
                    QTAILQ_INSERT_TAIL(&lm->line_resource[channel][lun].victim_line_list, line, victim_entry[channel][lun]);
                    line->inserted_to_victim_queue[channel][lun] = true;
                    lm->line_resource[channel][lun].victim_line_cnt++;
                }
            }else{
                // write_log("debug 8.3.8.0\n");
            }
        }
    }
}

static void mark_page_valid(struct ssd *ssd, struct ppa *ppa)
{
    struct nand_block *blk = NULL;
    struct nand_page *pg = NULL;
    struct line *line;

    /* update page status */
    pg = get_pg(ssd, ppa);
    if (pg->status != PG_FREE){
        write_log("ppa = %"PRIu64", lpa = %"PRIu64", pg->status = %d", ppa->ppa, get_rmap_ent(ssd, ppa), pg->status);
        
    }
    ftl_assert(pg->status == PG_FREE);
    pg->status = PG_VALID;

    /* update corresponding block status */
    blk = get_blk(ssd, ppa);
    ftl_assert(blk->vpc >= 0 && blk->vpc < ssd->sp.pgs_per_blk);
    blk->vpc++;

    /* update corresponding line status */
    line = get_line(ssd, ppa);
    ftl_assert(line->vpc >= 0 && line->vpc < line->pgs_per_line);
    line->vpc++;
}

static void mark_block_free(struct ssd *ssd, struct ppa *ppa)
{
    struct ssdparams *spp = &ssd->sp;
    struct nand_block *blk = get_blk(ssd, ppa);
    struct nand_page *pg = NULL;
    struct block_num *block_num = NULL;
    for (int i = 0; i < spp->pgs_per_blk; i++) {
        /* reset page status */
        pg = &blk->pg[i];
        ftl_assert(pg->nsecs == spp->secs_per_pg);
        pg->status = PG_FREE;
    }

    block_num = g_malloc0(sizeof(struct block_num));
    block_num->block_num = ppa->g.blk;
    // write_log("Breaking i: %d, j: %d, channel: %d, lun: %d, block: %d\n", i, j, channel, lun, block_num->block_num);
    // fflush(NULL);
    QTAILQ_INSERT_TAIL(&(ssd->channel_mgmt.channel[ppa->g.ch].lun[ppa->g.lun].free_block_list), block_num, entry);
    ssd->channel_mgmt.channel[ppa->g.ch].lun[ppa->g.lun].free_blocks_cnt++;

    /* reset block status */
    ftl_assert(blk->npgs == spp->pgs_per_blk);
    blk->ipc = 0;
    blk->vpc = 0;
    blk->erase_cnt++;
    blk->status = USE_FREE;
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
    struct ssdparams *spp = &ssd->sp;
    // write_log("debug gcwp 1\n");
    
    uint64_t lpn = get_rmap_ent(ssd, old_ppa);
    // write_log("debug gcwp 2\n");
    
    ftl_assert(valid_lpn(ssd, lpn));
    // For garbage collection, we just assign stream 0
    // since our previous guess of lifetime failed
    // write_log("debug gcwp 3\n");
    
    new_ppa = get_new_page(ssd, spp->gc_stream_id);
    /* update maptbl */
    // write_log("debug gcwp 4\n");
    
    set_maptbl_ent(ssd, lpn, &new_ppa);
    // write_log("debug gcwp 5\n");
    
    /* update rmap */
    set_rmap_ent(ssd, lpn, &new_ppa);
    // write_log("debug gcwp 6\n");
    
    mark_page_valid(ssd, &new_ppa);
    // write_log("debug gcwp 7\n");
    
    /* need to advance the write pointer here */
    ssd_advance_write_pointer(ssd, spp->gc_stream_id, lpn);
    // write_log("debug gcwp 8\n");
    
    if (ssd->sp.enable_gc_delay) {
        struct nand_cmd gcw;
        gcw.type = GC_IO;
        gcw.cmd = NAND_WRITE;
        gcw.stime = 0;
        ssd_advance_status(ssd, &new_ppa, &gcw);
    }
    // write_log("debug gcwp 9\n");
    
    /* advance per-ch gc_endtime as well */
#if 0
    new_ch = get_ch(ssd, &new_ppa);
    new_ch->gc_endtime = new_ch->next_ch_avail_time;
#endif
    // write_log("debug gcwp 10\n");
    
    new_lun = get_lun(ssd, &new_ppa);
    new_lun->gc_endtime = new_lun->next_lun_avail_time;
    // write_log("debug gcwp 11\n");
    
    return 0;
}

static struct line *find_smallest_vpc(struct ssd *ssd, int channel, int lun){
    struct line_mgmt *lm = &ssd->lm;
    struct line *line = NULL;
    struct line *temp_line = NULL;
    double min_vpc_ratio = 1.1;
    double temp_vpc_ratio;
    // write_log("finding smallest vpc for channel %d, lun %d\n", channel, lun);
    
    if (QTAILQ_EMPTY(&lm->line_resource[channel][lun].victim_line_list)) {
        return NULL;
    }
    QTAILQ_FOREACH(temp_line, &lm->line_resource[channel][lun].victim_line_list, victim_entry[channel][lun]){
        if (!temp_line->valid){
            write_log("bug line_id: %d\n", temp_line->id);
            fflush(NULL);
            assert(temp_line->valid);
        }
        temp_vpc_ratio = (double)temp_line->vpc / (double)temp_line->pgs_per_line;
        if (temp_vpc_ratio < min_vpc_ratio){
            line = temp_line;
            min_vpc_ratio = temp_vpc_ratio;
        }
    }
    // write_log("finished finding smallest vpc line %d for channel %d, lun %d\n", line->id, channel, lun);
    
    return line;
}

static struct line *select_victim_line(struct ssd *ssd, int channel, int lun, bool force)
{
    struct line_mgmt *lm = &ssd->lm;
    struct line *victim_line = NULL;
    int i, j;
    
    victim_line = find_smallest_vpc(ssd, channel, lun);
    // write_log("Finished finding smallest vpc\n");
    
    if (!victim_line) {
        return NULL;
    }
    // write_log("gc channel: %d, lun: %d, min_vpc: %d\n", channel, lun, victim_line->vpc);
    if (!force && victim_line->ipc < victim_line->pgs_per_line / 8) {
        return NULL;
    }

    for (i = 0; i < victim_line->total_channels; i++){
        channel = victim_line->channel_list[i];
        for (j = 0; j < victim_line->total_luns[i]; j++){
            lun = victim_line->lun_list[i][j];
            if (victim_line->inserted_to_victim_queue[channel][lun]){
                QTAILQ_REMOVE(&lm->line_resource[channel][lun].victim_line_list, victim_line, victim_entry[channel][lun]);
                lm->line_resource[channel][lun].victim_line_cnt--;
                victim_line->inserted_to_victim_queue[channel][lun] = false;
            }
        }
    }
    victim_line->valid = false;
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

static void mark_line_free(struct ssd *ssd, struct line *line)
{
    int i;
    struct ssdparams *spp = &ssd->sp;
    line->ipc = 0;
    line->vpc = 0;
    line->use = USE_FREE;

    for (i = 0; i < spp->nchs; i++){
        g_free(line->full_entry[i]);
        g_free(line->victim_entry[i]);
        g_free(line->inserted_to_full_queue[i]);
        g_free(line->inserted_to_victim_queue[i]);
    }
    
    for (i = 0; i < line->total_channels; i++){
        g_free(line->lun_list[i]);
        g_free(line->block_list[i]);
    }

    g_free(line->full_entry);
    g_free(line->victim_entry);
    g_free(line->inserted_to_full_queue);
    g_free(line->inserted_to_victim_queue);
    g_free(line->channel_list);
    g_free(line->lun_list);
    g_free(line->block_list);
    g_free(line->total_luns);

    g_free(line);
}

static struct ppa select_victim_block(struct ssd *ssd, int channel, int lun){
    // Choose a victim block for the given channel and lun.

    int i;
    struct ssdparams *spp = &ssd->sp;
    struct ppa ppa;
    int victim_id = -1;
    int victim_valid_count = spp->pgs_per_blk + 1;
    struct nand_block *blk;
    ppa.g.ch = channel;
    ppa.g.lun = lun;
    ppa.g.pl = 0;
    for (i = 0; i < spp->blks_per_lun; i++){
        // Choose the block with smallest number of pages
        ppa.g.blk = i;
        blk = get_blk(ssd, &ppa);
        if (blk->status == USE_FULL && blk->vpc < victim_valid_count){
            victim_valid_count = blk->vpc;
            victim_id = i;
        }
    }
    ppa.g.blk = victim_id;
    if (victim_id == -1){
        ppa.ppa = INVALID_PPA;
    }
    return ppa;
}

/*
static struct ppa select_overall_victim_block(struct ssd *ssd){
    int i, j, blockid;
    struct ssdparams *spp = &ssd->sp;
    struct ppa ppa;
    struct ppa result;
    result.ppa = INVALID_PPA;
    ppa.g.pl = 0;
    int smallest_count = spp->pgs_per_blk + 1;
    for (i = 0; i < spp->nchs; i++){
        for (j = 0; j < spp->luns_per_ch; j++){
            ppa = select_victim_block(ssd, i, j);
            if (ppa.ppa == INVALID_PPA){
                continue;
            }
            if (get_blk(ssd, &ppa)->vpc < smallest_count){
                result = ppa;
                smallest_count = get_blk(ssd, &ppa)->vpc;
            }
        }
    }
    return result;
}
*/
static int do_gc(struct ssd *ssd, bool force)
{
    struct line *victim_line = NULL;
    struct ssdparams *spp = &ssd->sp;
    struct nand_lun *lunp;
    struct ppa ppa;
    int ch, lun;
    int loopi, loopj;
    int i, j;
    int result = -1;
    int channel_end = spp->gc_start_channel + spp->nchs;
    int lun_end = spp->gc_start_lun + spp->luns_per_ch;
    ppa.ppa = INVALID_PPA;
    
    if (spp->enable_partial_gc){
        if (!should_gc_channel_lun(ssd, spp->gc_start_channel, spp->gc_start_lun, force)){
            spp->gc_start_lun++;
            if (spp->gc_start_lun == spp->luns_per_ch){
                spp->gc_start_lun = 0;
                spp->gc_start_channel = (spp->gc_start_channel + 1) % spp->nchs;
            }
            return result;
        }
        ppa = select_victim_block(ssd, spp->gc_start_channel, spp->gc_start_lun);
        assert(ppa.ppa != INVALID_PPA);
        write_log("[8, %d, %d, %d, %d, %d, %d, %d]\n", get_line(ssd, &ppa)->id, get_blk(ssd, &ppa)->ipc, get_blk(ssd, &ppa)->vpc, ppa.g.ch, ppa.g.lun, ppa.g.blk, force);
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
        spp->gc_start_lun++;
        if (spp->gc_start_lun == spp->luns_per_ch){
            spp->gc_start_lun = 0;
            spp->gc_start_channel = (spp->gc_start_channel + 1) % spp->nchs;
        }
        return result;
    }

    for (loopi = spp->gc_start_channel; loopi < channel_end; loopi++){
        i = loopi % spp->nchs;
        for (loopj = spp->gc_start_lun; loopj < lun_end; loopj++){
            j = loopj % spp->luns_per_ch;
            if (!should_gc_channel_lun(ssd, i, j, force)){
                continue;
            }
            if (spp->enable_partial_gc){
                
            }else{
                victim_line = select_victim_line(ssd, i, j, force);
                if (!victim_line) {
                    spp->gc_start_lun = (j + 1) % spp->luns_per_ch;
                    continue;
                }else{
                    result = 0;
                }

                // Line ID, Invalid count, Valid count, Total channels, Total luns, Victim line stream ID, Force
                write_log("[3, %d, %d, %d, %d, %d, %d, %d]\n", victim_line->id, victim_line->ipc, victim_line->vpc, victim_line->total_channels, victim_line->total_luns[i], victim_line->stream, force);

                /* copy back valid data */
                for (ch = 0; ch < victim_line->total_channels; ch++) {
                    ppa.g.ch = victim_line->channel_list[ch];
                    for (lun = 0; lun < victim_line->total_luns[ch]; lun++) {
                        ppa.g.lun = victim_line->lun_list[ch][lun]; // This should be put to ppa.g.ch and ppa.g.blk when we exploit chip level parallelism
                        ppa.g.pl = 0;
                        ppa.g.blk = victim_line->block_list[ch][lun];
                        // ftl_debug("GC-ing line:%d,ch:%d,lun:%d,blk:%d,ipc=%d,victim=%d,full=%d,free=%d,hostpgs=%"PRIu64",gcpgs=%"PRIu64"\n", victim_line->id, ch, lun, ppa.g.blk,
                        //   victim_line->ipc, ssd->lm.victim_line_cnt, ssd->lm.full_line_cnt,
                        //   ssd->lm.free_line_cnt, ssd->pages_from_host, ssd->pages_from_gc);
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
                assert(ppa.ppa != INVALID_PPA);
                /* update line status */
                // There might be a possibility that that block/channel/lun combination is used
                // So we pass the line (instead of a PPA)
                mark_line_free(ssd, victim_line);
                spp->gc_start_lun = (j + 1) % spp->luns_per_ch;
                if (j + 1 > spp->luns_per_ch){
                    spp->gc_start_channel = (i + 1) % spp->nchs;
                }
                return result;
            }
        }
    }

    return result;
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
    write_log("[6, %"PRIu64", %"PRIu64"]\n", start_lpn, end_lpn);
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
    // write_log("debug 1\n");
    
    NvmeRwCmd *rw = (NvmeRwCmd *) &req->cmd;
    NvmeNamespace *ns = req->ns;
    uint16_t control = le16_to_cpu(rw->control);
    uint32_t dsmgmt = le32_to_cpu(rw->dsmgmt);
    bool stream = control & NVME_RW_DTYPE_STREAMS;
    uint64_t page_death_time = 0;
    uint64_t passed_epoch_since_start = get_passed_epoch_since_start(ssd);
    // uint64_t uptime = get_uptime(ssd);
    struct stream_info *si;
    // struct stream_info *cmp_si;
    int offset = 0;
    uint16_t dspec = (dsmgmt >> 16) & 0xFFFF; //Stream ID
    // write_log("debug 2\n");
    

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
    // write_log("debug 3\n");
    
    uint64_t lba = req->slba;
    struct ssdparams *spp = &ssd->sp;
    int len = req->nlb;
    uint64_t start_lpn = lba / spp->secs_per_pg;
    uint64_t end_lpn = (lba + len - 1) / spp->secs_per_pg;
    uint64_t pages_written = (end_lpn - start_lpn) + 1;
    // Add WA tracker
    ssd->pages_from_host += pages_written;
    struct ppa ppa;
    uint64_t lpn;
    uint64_t curlat = 0, maxlat = 0;
    uint64_t chunk;
    uint64_t prediction;
    // write_log("debug 4\n");
    
    // uint64_t prev_starter = start_lpn;
    int stream_choice = 0;
    int i;
    int r;
    uint64_t ts = qemu_clock_get_us(QEMU_CLOCK_REALTIME);
    // double stream_min_lifetime;
    // double cur_stream_max_lifetime;
    //write_log("++++Write start LPN: %"PRIu64", end LPN: %"PRIu64", given stream: %d, now start++++\n", start_lpn, end_lpn, dspec);

    ////write_log("%s, opcode:%#x, start_sec:%#lx, size:%#lx, streamenabled:%d, dspec:%#x\n", __func__, rw->opcode, start_lpn * ssd->sp.secs_per_pg, (end_lpn - start_lpn + 1) * ssd->sp.secsz * ssd->sp.secs_per_pg, stream, dspec);
    // write_log("debug 5\n");
    
    if (end_lpn >= spp->tt_pgs) {
        ftl_err("start_lpn=%"PRIu64",tt_pgs=%d\n", start_lpn, ssd->sp.tt_pgs);
    }
    // write_log("debug 6\n");
    
    while (should_gc_high(ssd)) {
        /* perform GC here until !should_gc(ssd) */
        r = do_gc(ssd, true);
        if (r == -1)
            break;
    }
    // write_log("debug 7\n");
    
    // DZ Start
    // This is a write. Update death time and average.
    if (spp->death_time_prediction){
        set_latest_access_time(n, ssd, start_lpn, end_lpn, WRITE_OP);
    }
    // write_log("debug 8\n");
    
    // DZ End

    if (dspec > 0){
        stream_choice = dspec;
    } else if (spp->death_time_prediction){
        stream_choice = 0;
        chunk = lpn_to_chunk(start_lpn, n->pages_per_chunk);
        if (ssd->death_time_list[chunk].last_access_op != INITIAL_OP && ssd->death_time_list[chunk].last_access_op != WRITE_ONCE_OP){
            prediction = ssd->death_time_list[chunk].death_time_avg;
            page_death_time = passed_epoch_since_start + prediction;
            // Stream should have exponential accepted range
            // e.g. stream 1: DT b/t 0~1
            //      stream 2: DT b/t 1~3
            //      stream 3: DT b/t 3~7
            for (i = 0; i < spp->msl / 2; i++){
                prediction >>= 1;
                stream_choice += 1;
                if (prediction == 0){
                    break;
                }
            }

        }else if(ssd->death_time_list[chunk].last_access_op == WRITE_ONCE_OP){
            stream_choice = 0; // default stream
            //write_log("Addr: %"PRIu64", Chunk: %"PRIu64", First time written, no DT info, assigned to stream 0.\n", lpn, chunk);
        }
    }
    if (spp->enable_hetero_sbsize){
        stream_choice *= 2;
        if (pages_written < 8){
            offset = 0;
        }else{
            offset = 1;
        }
        stream_choice += offset;
    }

    si = &ssd->stream_info[stream_choice];

    for (lpn = start_lpn; lpn <= end_lpn; lpn++) {
        // write_log("debug 8.1\n");
        
        ppa = get_maptbl_ent(ssd, lpn);
        // write_log("debug 8.2\n");
        
        if (mapped_ppa(&ppa)) {
            /* update old page information first */
            // write_log("debug 8.3\n");
            
            mark_page_invalid(ssd, &ppa);
            // write_log("debug 8.4\n");
            
            set_rmap_ent(ssd, INVALID_LPN, &ppa);
        }
        // write_log("debug 8.5\n");
        
        // Update the write pointer earliest/latest death time for the target block
        if (page_death_time < si->earliest_death_time){
            si->earliest_death_time = page_death_time;
        }
        if (page_death_time > si->latest_death_time){
            si->latest_death_time = page_death_time;
        }
        // write_log("debug 10\n");
        
        /* new write */
        ppa = get_new_page(ssd, stream_choice);
        // write_log("debug 11\n");
        
        /* update maptbl */
        set_maptbl_ent(ssd, lpn, &ppa);
        // write_log("debug 12\n");
        
        /* update rmap */
        set_rmap_ent(ssd, lpn, &ppa);
        // write_log("debug 13\n");
        
        mark_page_valid(ssd, &ppa);
        // write_log("debug 14\n");
        
        /* need to advance the write pointer here */
        ssd_advance_write_pointer(ssd, stream_choice, lpn);
        // write_log("debug 15\n");
        
        struct nand_cmd swr;
        swr.type = USER_IO;
        swr.cmd = NAND_WRITE;
        swr.stime = req->stime;
        /* get latency statistics */
        curlat = ssd_advance_status(ssd, &ppa, &swr);
        maxlat = (curlat > maxlat) ? curlat : maxlat;
        // write_log("debug 16\n");
        
    }
    write_log("[1, %"PRIu64", %"PRIu64", %"PRIu64", %d]\n", start_lpn, end_lpn, ts, stream_choice);
    si->next_avail_time = maxlat;
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
            write_log("[2, %"PRIu64", %"PRIu64"]\n", start_lpn, end_lpn);
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

