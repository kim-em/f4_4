#!/usr/bin/env python3
"""Exact reconstruction of the whole C(F4,4) F-table in the normal-form gauge, over
k = Q(zeta_13)^+, written in the binary format of ../f4/fsym_bin.py (real basis of Q(zeta_52)^+,
so that ../f4/fsym_bin.Reader, fsym_k.py and the ../f4 verifiers read it unchanged).

Input: residues of the normal-form table (direct_f4.py) DIR/F_<p>_<w>.npy, for the primes P at the
roots w = 1,3,5,7,9,11 (one for each embedding of k) and for the held-out primes Q at all 12 real roots;
the primes in DIR/P.txt and DIR/Q.txt (those used for the archive: data/P.txt, data/Q.txt).  The block
layout (keys, row and column labels) is read from an existing table of the same layout,
f4_canon.PARTS ($F4_ARCHIVE/F4_4_F_part{1,2}.bin); only its index is used, not its entries.
Output: OUTDIR/F4_4_F_part1.bin, OUTDIR/F4_4_F_part2.bin.
Per block-aligned chunk: coefficients in the basis eta_1..eta_6 modulo p (6 x 6 solve), Garner CRT
(vectorised), rational reconstruction, check of every coefficient against the primes Q, check that
the residues of Q at w and at 25w agree (sigma_27-invariance, i.e. the entries lie in k), exact
conversion eta_j -> real basis (integer matrix), encoding with fsym_bin.encode_block.

usage: f4_full_recon.py DIR OUTDIR [--nproc n] [--target entries]
"""
import json, math, os, pickle, struct, sys, time
import multiprocessing as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "f4"))
import fsym_bin                       # read only
from lin import Fp, primitive_root_of_unity
from recon import ratrec, inv_mat
import cyc
import f4_canon as fc

D = sys.argv[1]
OUT = sys.argv[2]
NPROC = int(sys.argv[sys.argv.index("--nproc") + 1]) if "--nproc" in sys.argv else 16
TARGET = int(sys.argv[sys.argv.index("--target") + 1]) if "--target" in sys.argv else 300000
P = [int(x) for x in open(os.path.join(D, "P.txt")).read().split()]
Q = [int(x) for x in open(os.path.join(D, "Q.txt")).read().split()]
KREP = [1, 3, 5, 7, 9, 11]
PARTNER = {1: 25, 3: 23, 5: 21, 7: 19, 9: 17, 11: 15}
ETA2REAL = [[int(x) for x in cyc.to_real(e)] for e in cyc.ETA_BASIS]     # 6 x 12 integers
where = fc.readers()
keys = sorted(where)
sizes = [len(where[k].index[k]["rows"]) * len(where[k].index[k]["cols"]) for k in keys]
offs = np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)
VINV = {}
for p in P + Q:
    z = primitive_root_of_unity(p, 52, 1)
    V = [[(pow(z, 4 * j * w, p) + pow(z, -4 * j * w, p)) % p for j in range(1, 7)] for w in KREP]
    VINV[p] = np.array(inv_mat(V, p), dtype=np.int64)
PROD = [1]
for p in P:
    PROD.append(PROD[-1] * p)
M = PROD[-1]


def arr(p, w, a, b):
    return np.load(os.path.join(D, f"F_{p}_{w}.npy"), mmap_mode="r")[a:b]


def coefficients(p, a, b):
    vals = np.array([arr(p, w, a, b) for w in KREP], dtype=np.int64)
    return Fp(p).mm(VINV[p], vals)                                     # 6 x n


def work(rng):
    i0, i1, part = rng
    a, b = int(offs[i0]), int(offs[i1])
    co = [coefficients(p, a, b) for p in P]
    dig = [co[0]]
    for i in range(1, len(P)):
        p = P[i]
        t = dig[i - 1].copy()
        for j in range(i - 2, -1, -1):
            t = (t * (P[j] % p) + dig[j]) % p
        inv = pow(PROD[i] % p, -1, p)
        dig.append((co[i] - t) % p * inv % p)
    chk = {q: coefficients(q, a, b) for q in Q}
    bad = 0
    tau_bad = 0
    for q in Q:
        for w, w2 in PARTNER.items():
            tau_bad += int((arr(q, w, a, b) != arr(q, w2, a, b)).sum())
    nz = np.zeros(dig[0].shape, dtype=bool)
    for dd in dig:
        nz |= dd != 0
    n = b - a
    eta = [[(0, 1)] * 6 for _ in range(n)]
    hmax = 0
    ks, ts = np.nonzero(nz)
    for k, t in zip(ks.tolist(), ts.tolist()):
        x = sum(int(dig[i][k, t]) * PROD[i] for i in range(len(P)))
        r = ratrec(x, M)
        if r is None:
            bad += 1; continue
        nu, de = r
        if any((nu - int(chk[q][k, t]) * de) % q for q in Q):
            bad += 1; continue
        row = list(eta[t]); row[k] = (nu, de); eta[t] = row
        hmax = max(hmax, abs(nu).bit_length(), de.bit_length())
    for q in Q:
        bad += int((chk[q][~nz] != 0).sum())
    # exact conversion to the real basis
    num = [[0] * 12 for _ in range(n)]
    den = [[1] * 12 for _ in range(n)]
    for t in range(n):
        row = eta[t]
        if not any(nu for nu, _ in row):
            continue
        L = 1
        for nu, de in row:
            if nu:
                L = L * de // math.gcd(L, de)
        ints = [nu * (L // de) for nu, de in row]
        for kk in range(12):
            s = sum(ints[j] * ETA2REAL[j][kk] for j in range(6))
            if s:
                g = math.gcd(s, L)
                num[t][kk], den[t][kk] = s // g, L // g
    path = os.path.join(OUT, f".chunk_{part}_{i0:06d}.tmp")
    ents = []
    with open(path, "wb") as f:
        for bi in range(i0, i1):
            key = keys[bi]
            ix = where[key].index[key]
            s, e = int(offs[bi]) - a, int(offs[bi + 1]) - a
            ent, buf = fsym_bin.encode_block(list(key), ix["rows"], ix["cols"], num[s:e], den[s:e])
            ent["_len"] = len(buf)
            ents.append(ent)
            f.write(buf)
    return part, i0, bad, tau_bad, hmax, int(nz.sum()), ents, path


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    part_of = lambda k: 1 if k[0] in ("0011", "0012", "1010") else 2
    ranges = []
    for part in (1, 2):
        idx = [i for i, k in enumerate(keys) if part_of(k) == part]
        # contiguous runs (keys sorted, so the part-1 blocks form contiguous runs by first leg)
        runs, s = [], idx[0]
        for u, v in zip(idx, idx[1:] + [None]):
            if v != u + 1:
                runs.append((s, u + 1)); s = v
        for (s, e) in runs:
            i0, acc = s, 0
            for i in range(s, e):
                acc += sizes[i]
                if acc >= TARGET:
                    ranges.append((i0, i + 1, part)); i0, acc = i + 1, 0
            if i0 < e:
                ranges.append((i0, e, part))
    print(f"{len(keys)} blocks, {int(offs[-1])} entries, {len(ranges)} chunks; primes {P}, held out {Q}; modulus {M.bit_length()} bits", flush=True)
    with mp.get_context("fork").Pool(NPROC) as pool:
        res = []
        for r in pool.imap_unordered(work, ranges, chunksize=1):
            res.append(r)
            if len(res) % 20 == 0:
                print(f"  {len(res)}/{len(ranges)} chunks ({time.time() - t0:.0f} s)", flush=True)
    bad = sum(r[2] for r in res); tau = sum(r[3] for r in res); hmax = max(r[4] for r in res); nnz = sum(r[5] for r in res)
    print(f"reconstruction failures {bad}; residues of the held-out primes differing between w and 25w: {tau}; "
          f"largest numerator/denominator (eta basis) {hmax} bits; nonzero eta coefficients {nnz} ({time.time() - t0:.0f} s)", flush=True)
    if bad or tau:
        print("NOT written"); return 1
    # The header text below is part of the archived files (kept verbatim so that the output is
    # identical byte for byte; the paths in it are those of the machine on which the archive was made).
    meta = dict(field="k = Q(zeta_13)^+ inside Q(zeta_52)^+, zeta = exp(2 pi i/52); every entry lies in k",
                basis="real12: a0 + sum_k a_k (zeta^k + zeta^-k)",
                convention=("export: phi^{ZW->T}_j o (phi^{XY->Z}_i (x) 1) = sum F[(X,Y,W,T)][(Z,i,j),(Z',k,l)] "
                            "phi^{XZ'->T}_l o (1 (x) phi^{YW->Z'}_k)"),
                gauge="normal form of fields/f4_canon.py (plan fields/data/F4_4_canon.plan.pkl); "
                      "equals g.F for g in fields/data-large/F4_4_gauge.txt and F the table of ../f4",
                category="C(F4,4)", objects=sorted({k[0] for k in keys}), unit="0000",
                primes=P, check=Q)
    for part in (1, 2):
        rs = sorted([r for r in res if r[0] == part], key=lambda r: r[1])
        entries, off = [], 0
        for r in rs:
            for ent in r[6]:
                L = ent.pop("_len")
                ent["offset"] = off
                off += L
                entries.append(ent)
        h = json.dumps(dict(meta, part=part, blocks=entries)).encode()
        path = os.path.join(OUT, f"F4_4_F_part{part}.bin")
        with open(path, "wb") as f:
            f.write(fsym_bin.MAGIC)
            f.write(struct.pack("<Q", len(h)))
            f.write(h)
            for r in rs:
                with open(r[7], "rb") as g:
                    while True:
                        buf = g.read(1 << 26)
                        if not buf:
                            break
                        f.write(buf)
                os.remove(r[7])
        print(f"wrote {path}: {len(entries)} blocks, {os.path.getsize(path) / 1e9:.2f} GB", flush=True)
    print(f"total {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
