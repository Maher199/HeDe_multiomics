#!/usr/bin/env python3
"""CpG island detector (Gardiner-Garden & Frommer 1987), per-chromosome.

Criteria over a sliding WIN=200 window: GC% >= MIN_GC and Obs/Exp CpG >= MIN_OE.
Islands = maximal runs of qualifying windows with total length >= MIN_LEN.

Memory model: one chromosome at a time. For a chromosome of length L we hold
three int64 cumsum arrays (~3*8*L bytes). For the largest GRCr8 chromosome
(~280 Mb) that is ~6.7 GB transiently; we free it before the next chromosome.
Run solo on a >=16 GB machine. Skips NW_/NT_ scaffolds.

Output BED: chrom, start(0-based), end, cpgCount, gcPct, oeRatio.
"""
import sys, gzip, gc
import numpy as np

FASTA = sys.argv[1] if len(sys.argv) > 1 else "/mnt/shared-workspace/shared/rn_GRCr8.fa"
OUT   = sys.argv[2] if len(sys.argv) > 2 else "/workspace/hede_followup/ref/GRCr8_cgi.bed"
# Standard CpG island definition (UCSC/EMBOSS cpgplot; Takai-Jones & Gardiner-Garden):
# sliding 200-bp window, GC>=50%, Obs/Exp CpG>=0.60; merged islands must be >=500 bp
# (the >=500 bp island threshold yields the canonical mammalian CGI set and excludes
#  degenerate single-window artifacts).
WIN, MIN_GC, MIN_OE, MIN_LEN = 200, 0.50, 0.60, 500
MIN_CG = 8          # minimum CpG dinucleotides per qualifying window (excludes C/G homopolymers)
MIN_BASE = 10       # minimum C and minimum G count per window (excludes degenerate windows)
C_, G_ = 67, 71

def iter_chroms(path):
    name, chunks = None, []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]; chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)

def islands_for_chrom(seq):
    n = len(seq)
    if n < WIN:
        return []
    a = np.frombuffer(seq.upper().encode("ascii"), dtype=np.uint8)
    # int32 cumsums are safe: max count per chromosome < 2^31
    c = (a == C_); g = (a == G_)
    cg = np.zeros(n, dtype=bool); cg[:-1] = c[:-1] & g[1:]
    cc = np.concatenate(([0], np.cumsum(c, dtype=np.int32)))
    cg_ = np.concatenate(([0], np.cumsum(g, dtype=np.int32)))
    ccpg = np.concatenate(([0], np.cumsum(cg, dtype=np.int32)))
    del a, c, g, cg; gc.collect()
    # compute qualifying windows in blocks to bound the per-window arrays
    nwin = n - WIN + 1
    qual = np.empty(nwin, dtype=bool)
    BLK = 5_000_000
    for s in range(0, nwin, BLK):
        e = min(s + BLK, nwin)
        st = np.arange(s, e, dtype=np.int64)
        c_cnt = (cc[st+WIN] - cc[st]).astype(np.int64)
        g_cnt = (cg_[st+WIN] - cg_[st]).astype(np.int64)
        cp_cnt = (ccpg[st+WIN] - ccpg[st]).astype(np.int64)
        gc_frac = (c_cnt + g_cnt) / WIN
        with np.errstate(divide="ignore", invalid="ignore"):
            oe = np.where((c_cnt > 0) & (g_cnt > 0), (cp_cnt * WIN) / (c_cnt * g_cnt), 0.0)
        qual[s:e] = (gc_frac >= MIN_GC) & (oe >= MIN_OE) & (cp_cnt >= MIN_CG) & \
                    (c_cnt >= MIN_BASE) & (g_cnt >= MIN_BASE)
    # maximal runs of qualifying windows
    dq = np.diff(qual.astype(np.int8))
    run_starts = list(np.where(dq == 1)[0] + 1)
    run_ends = list(np.where(dq == -1)[0])  # inclusive index of last True
    if qual[0]:
        run_starts = [0] + run_starts
    if qual[-1]:
        run_ends = run_ends + [len(qual) - 1]
    out = []
    for rs, re in zip(run_starts, run_ends):
        istart = int(rs)
        iend = int(re) + WIN  # exclusive
        ilen = iend - istart
        if ilen < MIN_LEN:
            continue
        c_c = int(cc[iend] - cc[istart]); g_c = int(cg_[iend] - cg_[istart])
        cp = int(ccpg[iend] - ccpg[istart])
        if c_c == 0 or g_c == 0:
            continue
        gcp = (c_c + g_c) / ilen
        oer = (cp * ilen) / (c_c * g_c)
        # apply Gardiner-Garden criteria to the final merged island span
        if gcp >= MIN_GC and oer >= MIN_OE and cp >= MIN_CG:
            out.append((istart, iend, cp, round(gcp, 4), round(oer, 4)))
    return out

def main():
    n_written = 0
    with open(OUT, "w") as out:
        for name, seq in iter_chroms(FASTA):
            if name.startswith(("NW_", "NT_")):
                del seq; gc.collect(); continue
            for (s, e, cp, gcp, oe) in islands_for_chrom(seq):
                out.write(f"{name}\t{s}\t{e}\t{cp}\t{gcp:.4f}\t{oe:.4f}\n")
                n_written += 1
            print(f"{name}\t{n_written}", flush=True)
            del seq; gc.collect()
    print(f"TOTAL islands: {n_written}", flush=True)

if __name__ == "__main__":
    main()
