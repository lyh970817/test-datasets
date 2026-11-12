# A stratum with a shared missing block

`genotypes/study_eur` from the parent directory with every call of seven samples set missing on chromosomes 1
and 2, to exercise the order in which variant and sample missingness thresholds are applied. A batch of
samples that failed on part of the array shares one block of missing calls: a single variant threshold
applied first removes the whole block to keep the samples, while the cumulative call-rate ladder removes the
samples at its first rung and keeps the variants.

| Path | Contents |
| --- | --- |
| `study_eur_missing_block.{bed,bim,fam}` | `study_eur` (3962 variants, 125 samples) with the calls of `EUR_S100`-`EUR_S106` set missing at the 615 chr1 and chr2 variants (312 + 303) |

The `.bim` and `.fam` are those of `study_eur`, byte for byte, and every other call of the `.bed` is unchanged.
The seven samples carry no planted defect in `study_eur`, so the block is their only one.

## Expected outcome

Measured with plink2 2.0.0a.6.9; the sample rates are over the whole fileset, and the two filter steps run over
the 3954 variants that pass `--maf 0.01`:

| Statistic | Value |
| --- | --- |
| variant F_MISS at the 615 block variants | 0.056-0.080 (7 of 125 samples, plus the random calls of the planted high-missingness samples) |
| sample F_MISS of `EUR_S100`-`EUR_S106` | 0.1562-0.1579 |
| variants left by `--geno 0.05` alone, then samples left by `--mind 0.05` | 3315 variants (the 615 block variants and 24 others removed), 123 samples (only `EUR_S040` and `EUR_S041` removed) |
| `study_eur` under the same two steps | 3930 variants, 123 samples |

On the ladder (`--geno` then `--mind` at 0.10, applied to the fileset each step leaves, then 0.09, ... 0.05),
the seven block samples fall at the first rung's sample step, F_MISS 0.156 > 0.10, and the block variants,
then missing in none of the remaining samples, survive every later rung.

## Reproduction

With Python 3 (standard library only) from the parent directory. A `.bed` is variant-major with `ceil(n / 4)`
bytes per variant after its 3-byte header, two bits per sample from the low bits up, `01` for a missing call:

```bash
python3 - genotypes/study_eur missing_block/study_eur_missing_block <<'PY'
import sys
src, out = sys.argv[1], sys.argv[2]
block_samples = {f"EUR_S{n:03d}" for n in range(100, 107)}
block_chroms = {"1", "2"}
fam = open(src + ".fam").read().splitlines()
bim = open(src + ".bim").read().splitlines()
n = len(fam)
idx = [i for i, line in enumerate(fam) if line.split()[1] in block_samples]
assert len(idx) == len(block_samples)
width = (n + 3) // 4
bed = bytearray(open(src + ".bed", "rb").read())
assert bed[:3] == b"\x6c\x1b\x01" and len(bed) == 3 + width * len(bim)
for v, line in enumerate(bim):
    if line.split("\t")[0] in block_chroms:
        for i in idx:
            pos = 3 + v * width + i // 4
            shift = (i % 4) * 2
            bed[pos] = (bed[pos] & ~(3 << shift) & 0xFF) | (1 << shift)
open(out + ".bed", "wb").write(bytes(bed))
open(out + ".bim", "w").write(open(src + ".bim").read())
open(out + ".fam", "w").write(open(src + ".fam").read())
PY
```
