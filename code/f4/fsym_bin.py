"""Compact binary format for large exact F-tables over the real subfield of Q(zeta_52).

An entry x is a0 + sum_{k=1}^{11} a_k (zeta^k + zeta^-k), zeta = exp(2 pi i/52), a_k rational.

File layout (one file per chunk, little endian):
  b"FSYMBIN1\\n"
  8 bytes: length h of the JSON header, then h bytes of UTF-8 JSON:
    {"field": ..., "basis": "real12", "convention": ..., "objects": [...], "unit": ...,
     "blocks": [{"key": [X,Y,W,T], "rows": [[Z,i,j],...], "cols": [[Z,k,l],...],
                 "offset": byte offset of the block data after the header, "nlimbs": n,
                 "dlimbs": m}, ...]}
  block data, for a block with r rows and c columns (e = r c entries, row-major):
    m uint32 limbs of the common denominator D of the block (least significant first),
    e x 12 int8 signs, then e x 12 x n uint32 limbs of the magnitudes |a_k D|.
So a_k = sign * magnitude / D exactly.  Limbs are 32-bit; magnitudes are padded to n limbs.

Reader (numpy only):  r = Reader(path);  r.keys();  r.block(key) -> (rows, cols, D, signs, limbs);
r.fractions(key) -> nested lists of 12-tuples of Fraction;  r.residues(key, p, cvals) -> the
matrix modulo p of x evaluated at a root (cvals = [1, c_1, ..., c_11] mod p).
"""
import json, struct
from fractions import Fraction
import numpy as np

MAGIC = b"FSYMBIN1\n"


def _limbs(n, nl):
    out = []
    for _ in range(nl):
        out.append(n & 0xFFFFFFFF)
        n >>= 32
    assert n == 0
    return out


def encode_block(key, rows, cols, num, den):
    """num, den: e x 12 nested lists of Python ints (den > 0).  Returns (index entry, bytes)."""
    from math import lcm
    e = len(rows) * len(cols)
    D = 1
    for dd in den:
        for d in dd:
            D = lcm(D, d)
    mags = [[abs(a) * (D // d) for a, d in zip(nr, dr)] for nr, dr in zip(num, den)]
    nl = max(1, (max((m.bit_length() for row in mags for m in row), default=0) + 31) // 32)
    dl = max(1, (D.bit_length() + 31) // 32)
    signs = np.array([[(a > 0) - (a < 0) for a in nr] for nr in num], dtype=np.int8).reshape(e, 12)
    L = np.zeros((e, 12, nl), dtype=np.uint32)
    for t, row in enumerate(mags):
        for k, m in enumerate(row):
            i = 0
            while m:
                L[t, k, i] = m & 0xFFFFFFFF
                m >>= 32
                i += 1
    buf = np.array(_limbs(D, dl), dtype=np.uint32).tobytes() + signs.tobytes() + L.tobytes()
    return dict(key=list(key), rows=[list(r) for r in rows], cols=[list(c) for c in cols], nlimbs=nl, dlimbs=dl), buf


def write_encoded(path, meta, encoded):
    """encoded: list of (index entry, bytes) in file order."""
    off = 0
    entries = []
    for ent, buf in encoded:
        ent = dict(ent, offset=off)
        off += len(buf)
        entries.append(ent)
    h = json.dumps(dict(meta, blocks=entries)).encode()
    with open(path, "wb") as f:
        f.write(MAGIC)
        f.write(struct.pack("<Q", len(h)))
        f.write(h)
        for _, buf in encoded:
            f.write(buf)


def write(path, meta, blocks):
    """blocks: list of (key, rows, cols, num, den)."""
    write_encoded(path, meta, [encode_block(*b) for b in blocks])


class Reader:
    def __init__(self, path):
        self.f = open(path, "rb")
        assert self.f.read(len(MAGIC)) == MAGIC
        (h,) = struct.unpack("<Q", self.f.read(8))
        self.meta = json.loads(self.f.read(h).decode())
        self.base = len(MAGIC) + 8 + h
        self.index = {tuple(b["key"]): b for b in self.meta["blocks"]}

    def keys(self):
        return list(self.index)

    def block(self, key):
        b = self.index[tuple(key)]
        e = len(b["rows"]) * len(b["cols"])
        nl, dl = b["nlimbs"], b["dlimbs"]
        self.f.seek(self.base + b["offset"])
        dlimbs = np.frombuffer(self.f.read(4 * dl), dtype=np.uint32)
        D = sum(int(x) << (32 * i) for i, x in enumerate(dlimbs))
        signs = np.frombuffer(self.f.read(e * 12), dtype=np.int8).reshape(e, 12)
        limbs = np.frombuffer(self.f.read(4 * e * 12 * nl), dtype=np.uint32).reshape(e, 12, nl)
        return [tuple(r) for r in b["rows"]], [tuple(c) for c in b["cols"]], D, signs, limbs

    def fractions(self, key):
        rows, cols, D, signs, limbs = self.block(key)
        nl = limbs.shape[2]
        vals = []
        for t in range(len(rows) * len(cols)):
            vals.append(tuple(Fraction(int(signs[t, k]) * sum(int(limbs[t, k, i]) << (32 * i) for i in range(nl)), D)
                              for k in range(12)))
        c = len(cols)
        return rows, cols, [vals[i * c:(i + 1) * c] for i in range(len(rows))]

    def residues(self, key, p, cvals):
        """Matrix of the block modulo p at a root, cvals[k] = value of c_k (cvals[0] = 1)."""
        rows, cols, D, signs, limbs = self.block(key)
        nl = limbs.shape[2]
        acc = np.zeros(limbs.shape[:2], dtype=np.int64)
        base = (1 << 32) % p
        for i in range(nl - 1, -1, -1):              # Horner in the limb base 2^32
            acc = (acc * base + (limbs[:, :, i].astype(np.int64) % p)) % p
        acc = acc * (signs.astype(np.int64) % p) % p
        assert D % p, "p divides a denominator"
        cv = np.array([int(x) % p for x in cvals], dtype=np.int64)
        v = np.zeros(acc.shape[0], dtype=np.int64)
        for k in range(12):
            v = (v + acc[:, k] * cv[k]) % p
        v = v * pow(D, -1, p) % p
        return np.array(v).reshape(len(rows), len(cols))
