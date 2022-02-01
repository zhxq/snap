#include "./nvme.h"

// Methods for multi-stream support is adapted from here:
// https://github.com/multi-stream/qemu/commit/9888fffbf3695f7de4170e97f49231ea18fa8cd4#

int nvme_check_sqid(FemuCtrl *n, uint16_t sqid)
{
    return sqid <= n->num_io_queues && n->sq[sqid] != NULL ? 0 : -1;
}

int nvme_check_cqid(FemuCtrl *n, uint16_t cqid)
{
    return cqid <= n->num_io_queues && n->cq[cqid] != NULL ? 0 : -1;
}

void nvme_inc_cq_tail(NvmeCQueue *cq)
{
    cq->tail++;
    if (cq->tail >= cq->size) {
        cq->tail = 0;
        cq->phase = !cq->phase;
    }
}

void nvme_inc_sq_head(NvmeSQueue *sq)
{
    sq->head = (sq->head + 1) % sq->size;
}

void nvme_update_sq_tail(NvmeSQueue *sq)
{
    if (sq->db_addr_hva) {
        sq->tail = *((uint32_t *)sq->db_addr_hva);
        return;
    }

    if (sq->db_addr) {
        nvme_addr_read(sq->ctrl, sq->db_addr, &sq->tail, sizeof(sq->tail));
    }
}

void nvme_update_cq_head(NvmeCQueue *cq)
{
    if (cq->db_addr_hva) {
        cq->head = *(uint32_t *)(cq->db_addr_hva);
        return;
    }

    if (cq->db_addr) {
        nvme_addr_read(cq->ctrl, cq->db_addr, &cq->head, sizeof(cq->head));
    }
}

uint8_t nvme_cq_full(NvmeCQueue *cq)
{
    nvme_update_cq_head(cq);

    return (cq->tail + 1) % cq->size == cq->head;
}

uint8_t nvme_sq_empty(NvmeSQueue *sq)
{
    return sq->head == sq->tail;
}

uint64_t *nvme_setup_discontig(FemuCtrl *n, uint64_t prp_addr, uint16_t
                               queue_depth, uint16_t entry_size)
{
    uint16_t prps_per_page = n->page_size >> 3;
    uint64_t prp[prps_per_page];
    uint16_t total_prps = DIV_ROUND_UP(queue_depth * entry_size, n->page_size);
    uint64_t *prp_list = g_malloc0(total_prps * sizeof(*prp_list));

    for (int i = 0; i < total_prps; i++) {
        if (i % prps_per_page == 0 && i < total_prps - 1) {
            if (!prp_addr || prp_addr & (n->page_size - 1)) {
                g_free(prp_list);
                return NULL;
            }
            nvme_addr_write(n, prp_addr, (uint8_t *)&prp, sizeof(prp));
            prp_addr = le64_to_cpu(prp[prps_per_page - 1]);
        }
        prp_list[i] = le64_to_cpu(prp[i % prps_per_page]);
        if (!prp_list[i] || prp_list[i] & (n->page_size - 1)) {
            g_free(prp_list);
            return NULL;
        }
    }

    return prp_list;
}

void nvme_set_error_page(FemuCtrl *n, uint16_t sqid, uint16_t cid, uint16_t
                         status, uint16_t location, uint64_t lba, uint32_t nsid)
{
    NvmeErrorLog *elp;

    elp = &n->elpes[n->elp_index];
    elp->error_count = n->error_count++;
    elp->sqid = sqid;
    elp->cid = cid;
    elp->status_field = status;
    elp->param_error_location = location;
    elp->lba = lba;
    elp->nsid = nsid;
    n->elp_index = (n->elp_index + 1) % n->elpe;
    ++n->num_errors;
}

uint16_t femu_nvme_rw_check_req(FemuCtrl *n, NvmeNamespace *ns, NvmeCmd *cmd,
                                NvmeRequest *req, uint64_t slba, uint64_t elba,
                                uint32_t nlb, uint16_t ctrl, uint64_t data_size,
                                uint64_t meta_size)
{

    if (elba > le64_to_cpu(ns->id_ns.nsze)) {
        nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_LBA_RANGE,
                            offsetof(NvmeRwCmd, nlb), elba, ns->id);
        return NVME_LBA_RANGE | NVME_DNR;
    }
    if (n->id_ctrl.mdts && data_size > n->page_size * (1 << n->id_ctrl.mdts)) {
        nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_INVALID_FIELD,
                            offsetof(NvmeRwCmd, nlb), nlb, ns->id);
        return NVME_INVALID_FIELD | NVME_DNR;
    }
    if (meta_size) {
        nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_INVALID_FIELD,
                            offsetof(NvmeRwCmd, control), ctrl, ns->id);
        return NVME_INVALID_FIELD | NVME_DNR;
    }
    if ((ctrl & NVME_RW_PRINFO_PRACT) && !(ns->id_ns.dps & DPS_TYPE_MASK)) {
        nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_INVALID_FIELD,
                            offsetof(NvmeRwCmd, control), ctrl, ns->id);
        /* Not contemplated in LightNVM for now */
        if (OCSSD(n)) {
            return 0;
        }
        return NVME_INVALID_FIELD | NVME_DNR;
    }
    if (!req->is_write && find_next_bit(ns->uncorrectable, elba, slba) < elba) {
        nvme_set_error_page(n, req->sq->sqid, cmd->cid, NVME_UNRECOVERED_READ,
                            offsetof(NvmeRwCmd, slba), elba, ns->id);
        return NVME_UNRECOVERED_READ;
    }

    return 0;
}

void nvme_free_sq(NvmeSQueue *sq, FemuCtrl *n)
{
    n->sq[sq->sqid] = NULL;
    g_free(sq->io_req);
    if (sq->prp_list) {
        g_free(sq->prp_list);
    }
    if (sq->sqid) {
        g_free(sq);
    }
}

uint16_t nvme_init_sq(NvmeSQueue *sq, FemuCtrl *n, uint64_t dma_addr, uint16_t
                      sqid, uint16_t cqid, uint16_t size, enum NvmeQueueFlags
                      prio, int contig)
{
    uint8_t stride = n->db_stride;
    int dbbuf_entry_sz = 1 << (2 + stride);
    AddressSpace *as = pci_get_address_space(&n->parent_obj);
    dma_addr_t sqsz = (dma_addr_t)size;
    NvmeCQueue *cq;

    sq->ctrl = n;
    sq->sqid = sqid;
    sq->size = size;
    sq->cqid = cqid;
    sq->head = sq->tail = 0;
    sq->phys_contig = contig;
    if (sq->phys_contig) {
        sq->dma_addr = dma_addr;
        sq->dma_addr_hva = (uint64_t)dma_memory_map(as, dma_addr, &sqsz, 0);
    } else {
        sq->prp_list = nvme_setup_discontig(n, dma_addr, size, n->sqe_size);
        if (!sq->prp_list) {
            return NVME_INVALID_FIELD | NVME_DNR;
        }
    }

    sq->io_req = g_malloc0(sq->size * sizeof(*sq->io_req));
    QTAILQ_INIT(&sq->req_list);
    QTAILQ_INIT(&sq->out_req_list);
    for (int i = 0; i < sq->size; i++) {
        sq->io_req[i].sq = sq;
        QTAILQ_INSERT_TAIL(&(sq->req_list), &sq->io_req[i], entry);
    }

    switch (prio) {
    case NVME_Q_PRIO_URGENT:
        sq->arb_burst = (1 << NVME_ARB_AB(n->features.arbitration));
        break;
    case NVME_Q_PRIO_HIGH:
        sq->arb_burst = NVME_ARB_HPW(n->features.arbitration) + 1;
        break;
    case NVME_Q_PRIO_NORMAL:
        sq->arb_burst = NVME_ARB_MPW(n->features.arbitration) + 1;
        break;
    case NVME_Q_PRIO_LOW:
    default:
        sq->arb_burst = NVME_ARB_LPW(n->features.arbitration) + 1;
        break;
    }

    if (sqid && n->dbs_addr && n->eis_addr) {
        sq->db_addr = n->dbs_addr + 2 * sqid * dbbuf_entry_sz;
        sq->db_addr_hva = n->dbs_addr_hva + 2 * sqid * dbbuf_entry_sz;
        sq->eventidx_addr = n->eis_addr + 2 * sqid * dbbuf_entry_sz;
        sq->eventidx_addr = n->eis_addr_hva + 2 * sqid + dbbuf_entry_sz;
        femu_debug("SQ[%d],db=%" PRIu64 ",ei=%" PRIu64 "\n", sqid, sq->db_addr,
                sq->eventidx_addr);
    }

    assert(n->cq[cqid]);
    cq = n->cq[cqid];
    QTAILQ_INSERT_TAIL(&(cq->sq_list), sq, entry);
    n->sq[sqid] = sq;

    return NVME_SUCCESS;
}

uint16_t nvme_init_cq(NvmeCQueue *cq, FemuCtrl *n, uint64_t dma_addr, uint16_t
                      cqid, uint16_t vector, uint16_t size, uint16_t
                      irq_enabled, int contig)
{
    cq->ctrl = n;
    cq->cqid = cqid;
    cq->size = size;
    cq->phase = 1;
    cq->irq_enabled = irq_enabled;
    cq->vector = vector;
    cq->head = cq->tail = 0;
    cq->phys_contig = contig;

    uint8_t stride = n->db_stride;
    int dbbuf_entry_sz = 1 << (2 + stride);
    AddressSpace *as = pci_get_address_space(&n->parent_obj);
    dma_addr_t cqsz = (dma_addr_t)size;

    if (cq->phys_contig) {
        cq->dma_addr = dma_addr;
        cq->dma_addr_hva = (uint64_t)dma_memory_map(as, dma_addr, &cqsz, 1);
    } else {
        cq->prp_list = nvme_setup_discontig(n, dma_addr, size,
                n->cqe_size);
        if (!cq->prp_list) {
            return NVME_INVALID_FIELD | NVME_DNR;
        }
    }

    QTAILQ_INIT(&cq->req_list);
    QTAILQ_INIT(&cq->sq_list);
    if (cqid && n->dbs_addr && n->eis_addr) {
        cq->db_addr = n->dbs_addr + (2 * cqid + 1) * dbbuf_entry_sz;
        cq->db_addr_hva = n->dbs_addr_hva + (2 * cqid + 1) * dbbuf_entry_sz;
        cq->eventidx_addr = n->eis_addr + (2 * cqid + 1) * dbbuf_entry_sz;
        cq->eventidx_addr_hva = n->eis_addr_hva + (2 * cqid + 1) * dbbuf_entry_sz;
        femu_debug("CQ, db_addr=%" PRIu64 ", eventidx_addr=%" PRIu64 "\n",
                cq->db_addr, cq->eventidx_addr);
    }
    msix_vector_use(&n->parent_obj, cq->vector);
    n->cq[cqid] = cq;

    return NVME_SUCCESS;
}

void nvme_free_cq(NvmeCQueue *cq, FemuCtrl *n)
{
    n->cq[cq->cqid] = NULL;
    msix_vector_unuse(&n->parent_obj, cq->vector);
    if (cq->prp_list) {
        g_free(cq->prp_list);
    }
    if (cq->cqid) {
        g_free(cq);
    }
}

void nvme_set_ctrl_name(FemuCtrl *n, const char *mn, const char *sn, int *dev_id)
{
    NvmeIdCtrl *id = &n->id_ctrl;
    char *subnqn;
    char serial[MN_MAX_LEN], dev_id_str[ID_MAX_LEN];

    memset(serial, 0, MN_MAX_LEN);
    memset(dev_id_str, 0, ID_MAX_LEN);
    strcat(serial, sn);

    sprintf(dev_id_str, "%d", *dev_id);
    strcat(serial, dev_id_str);
    (*dev_id)++;
    strpadcpy((char *)id->mn, sizeof(id->mn), mn, ' ');

    memset(n->devname, 0, MN_MAX_LEN);
    g_strlcpy(n->devname, serial, sizeof(serial));

    strpadcpy((char *)id->sn, sizeof(id->sn), serial, ' ');
    strpadcpy((char *)id->fr, sizeof(id->fr), "1.3", ' ');

    subnqn = g_strdup_printf("nqn.2021-05.org.femu:%s", serial);
    strpadcpy((char *)id->subnqn, sizeof(id->subnqn), subnqn, '\0');
}

uint16_t nvme_dir_receive(FemuCtrl *n, NvmeCmd *cmd, NvmeCqe *cqe)
{
    NvmeRwCmd *rw = (NvmeRwCmd *)cmd;
    uint32_t nsid = le32_to_cpu(cmd->nsid);
    uint64_t prp1 = le64_to_cpu(rw->prp1);
    uint64_t prp2 = le64_to_cpu(rw->prp2);
    uint32_t numd  = le32_to_cpu(cmd->cdw10);
    uint32_t dw11  = le32_to_cpu(cmd->cdw11);
    uint32_t dw12  = le32_to_cpu(cmd->cdw12);
    //uint16_t dspec = (dw11 >> 16) & 0xFFFF;
    uint8_t  dtype  = (dw11 >> 8) & 0xFF;
    uint8_t  doper  = dw11 & 0xFF;
    uint16_t nsr;
    uint32_t result = 0;
    NvmeNamespace *ns;

    if (nsid != 0xffffffff && (nsid == 0 || nsid > n->num_namespaces)) {
        return NVME_INVALID_NSID | NVME_DNR;
    }
    if (nsid == 0xffffffff)
        ns = &n->namespaces[0];
    else
        ns = &n->namespaces[nsid - 1];

    if (!(n->id_ctrl.oacs & NVME_OACS_DIR))
        return NVME_INVALID_OPCODE;

    switch (dtype) {
    case NVME_DIR_TYPE_IDENTIFY:
        switch (doper) {
        case NVME_DIR_RCV_ID_OP_PARAM:
            if (((numd + 1) * 4) < sizeof(*(ns->id_dir)))
                return dma_read_prp(n,
                       (uint8_t *)ns->id_dir,
                       (numd + 1) * 4,
                       prp1, prp2);
            else
                return dma_read_prp(n,
                       (uint8_t *)ns->id_dir,
                       sizeof(*(ns->id_dir)),
                       prp1, prp2);
            break;
        default:
            return NVME_INVALID_FIELD;
        }
        break;
    case NVME_DIR_TYPE_STREAMS:
        switch (doper) {
        case NVME_DIR_RCV_ST_OP_PARAM:
            ns->str_ns_param->msl = n->str_sys_param->msl;
            ns->str_ns_param->nssa = n->str_sys_param->nssa;
            ns->str_ns_param->nsso = n->str_sys_param->nsso;
            if (((numd + 1) * 4) < sizeof(*(ns->str_ns_param)))
                return dma_read_prp(n,
                       (uint8_t *)ns->str_ns_param,
                       (numd + 1) * 4,
                       prp1, prp2);
            else
                return dma_read_prp(n,
                       (uint8_t *)ns->str_ns_param,
                       sizeof(*(ns->str_ns_param)),
                       prp1, prp2);
            break;
        case NVME_DIR_RCV_ST_OP_STATUS:
            if (((numd + 1) * 4) < sizeof(*(ns->str_ns_stat)))
                return dma_read_prp(n,
                       (uint8_t *)ns->str_ns_stat,
                       (numd + 1) * 4,
                       prp1, prp2);
            else
                return dma_read_prp(n,
                       (uint8_t *)ns->str_ns_stat,
                       sizeof(*(ns->str_ns_stat)),
                       prp1, prp2);
            break;
        case NVME_DIR_RCV_ST_OP_RESOURCE:
            if (ns->str_ns_param->nsa)
                return NVME_INVALID_FIELD;
            nsr = dw12 & 0xFFFF;
            if (nsr > n->str_sys_param->nssa)
                nsr = n->str_sys_param->nssa;
            ns->str_ns_param->nsa = nsr;
            n->str_sys_param->nssa -= nsr;
            result = cpu_to_le32(nsr);
            break;
        default:
            return NVME_INVALID_FIELD;
        }
        break;
    default:
        return NVME_INVALID_FIELD;
    }

    cqe->n.result = result;
    return NVME_SUCCESS;
}

int nvme_found_in_str_list(NvmeDirStrNsStat *str_ns_stat, uint16_t dspec)
{
    int i;

    if (str_ns_stat->cnt == 0)
        return -1;

    for (i=0; i<str_ns_stat->cnt; i++)
        if (str_ns_stat->id[i] == dspec)
            return i;
    return -1;

}

void nvme_add_to_str_list(NvmeDirStrNsStat *str_ns_stat, uint16_t dspec)
{
   str_ns_stat->id[str_ns_stat->cnt] = dspec;
   str_ns_stat->cnt++;
}

void nvme_del_from_str_list(NvmeDirStrNsStat *str_ns_stat, int pos)
{
   int i;

   if (str_ns_stat->cnt == 0)
        return;
   str_ns_stat->cnt--;
   for (i=pos; i<str_ns_stat->cnt; i++)
        str_ns_stat->id[i] = str_ns_stat->id[i+1];
   str_ns_stat->id[str_ns_stat->cnt] = 0;
}

void nvme_update_str_stat(FemuCtrl *n, NvmeNamespace *ns, uint16_t dspec)
{
    NvmeDirStrNsStat *st = ns->str_ns_stat;
    if (dspec == 0) /* skip if normal write */
        return;
    if (nvme_found_in_str_list(st, dspec) < 0) { /* not found */
        /* delete the first if max out */
        if (n->str_sys_param->nsso == n->str_sys_param->msl) {
        ns->str_ns_param->nso--;
        n->str_sys_param->nsso--;
        nvme_del_from_str_list(st, 0);
        }
        nvme_add_to_str_list(st, dspec);
        ns->str_ns_param->nso++;
        n->str_sys_param->nsso++;
    }
   return;
}

uint16_t nvme_dir_send(FemuCtrl *n, NvmeCmd *cmd)
{
    int i;
    uint32_t nsid = le32_to_cpu(cmd->nsid);
    //uint64_t prp1 = le64_to_cpu(cmd->prp1);
    //uint64_t prp2 = le64_to_cpu(cmd->prp2);
    //uint32_t numd = le32_to_cpu(cmd->cdw10);
    uint32_t dw11  = le32_to_cpu(cmd->cdw11);
    uint32_t dw12  = le32_to_cpu(cmd->cdw12);
    uint16_t dspec = (dw11 >> 16) & 0xFFFF;
    uint8_t  dtype  = (dw11 >> 8) & 0xFF;
    uint8_t  doper  = dw11 & 0xFF;
    uint8_t  tdtype;
    uint8_t  endir;
    NvmeNamespace *ns;

    if (nsid != 0xffffffff && (nsid == 0 || nsid > n->num_namespaces)) {
        return NVME_INVALID_NSID | NVME_DNR;
    }
    if (nsid == 0xffffffff)
        ns = &n->namespaces[0];
    else
        ns = &n->namespaces[nsid - 1];

    if (!(n->id_ctrl.oacs & NVME_OACS_DIR))
        return NVME_INVALID_OPCODE;

    switch (dtype) {
    case NVME_DIR_TYPE_IDENTIFY:
        switch (doper) {
        case NVME_DIR_SND_ID_OP_ENABLE:
            tdtype = (dw12 >> 8) & 0xFF;
            endir = dw12 & NVME_DIR_ENDIR;
            if (tdtype == NVME_DIR_TYPE_STREAMS) {
                if (endir)
                    ns->id_dir->dir_enable[0] |= NVME_DIR_IDF_STREAMS;
                else
                    ns->id_dir->dir_enable[0] &= ~NVME_DIR_IDF_STREAMS;
            }
            break;
        default:
            return NVME_INVALID_FIELD;
        }
        break;
    case NVME_DIR_TYPE_STREAMS:
        switch (doper) {
        case NVME_DIR_SND_ST_OP_REL_ID:
            if ((i = nvme_found_in_str_list(ns->str_ns_stat, dspec)) >= 0) {
                ns->str_ns_param->nso--;
                n->str_sys_param->nsso--;
                nvme_del_from_str_list(ns->str_ns_stat, i);
            }
            break;
        case NVME_DIR_SND_ST_OP_REL_RSC:
            ns->str_ns_param->nsa = 0;
            n->str_sys_param->nssa = n->str_sys_param->msl;
            break;
        default:
            return NVME_INVALID_FIELD;
        }
        break;
    default:
        return NVME_INVALID_FIELD;
    }

    return NVME_SUCCESS;
}