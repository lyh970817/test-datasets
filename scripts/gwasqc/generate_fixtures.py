#!/usr/bin/env python3
"""Deterministic generator for the gwasqc genotype-QC test fixtures.

Regenerates every file under ``data/gwasqc/`` from scratch, including the
README that documents the planted defects.  Pure Python standard library (no
numpy / no R), single fixed seed, so repeated runs are byte-identical.

Usage (from the repository root of this branch)::

    python3 scripts/gwasqc/generate_fixtures.py

Options::

    --out DIR     output directory (default: data/gwasqc)
    --seed INT    master seed (default: 20260920)

Design notes
------------
* Allele frequencies follow a two-level Balding-Nichols model: an ancestral
  frequency per variant, differentiated across five super-populations with
  Fst = 0.12, then across two sub-populations inside each super-population
  with Fst = 0.01.  Genotypes are binomial draws from the relevant frequency,
  so every variant is HWE-conforming unless a defect was planted into it.
* The reference panel carries autosomes only (that is all the ancestry PCA
  needs).  Study strata carry chr1-22, chrX (PAR1 / non-PAR / PAR2) and chrY.
* Study strata share variant IDs, positions and alleles with the reference on
  the autosomes, except for a documented ~1.5% strand-flipped subset whose
  study alleles are the complement of the reference alleles.
* Defects are planted into the study strata only; the reference panel is clean.
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

# --------------------------------------------------------------------------- #
# Parameters
# --------------------------------------------------------------------------- #

DEFAULT_SEED = 20260920

# GRCh38 primary assembly lengths (chr1-22).
CHROM_LEN_B38 = {
    1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259,
    6: 170805979, 7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422,
    11: 135086622, 12: 133275309, 13: 114364328, 14: 107043718, 15: 101991189,
    16: 90338345, 17: 83257441, 18: 80373285, 19: 58617616, 20: 64444167,
    21: 46709983, 22: 50818468,
}
AUTOSOME_TARGET = 3600
AUTOSOME_MIN_PER_CHROM = 60
VARIANT_MARGIN = 1_000_000  # keep variants away from telomeres

# GRCh38 pseudoautosomal boundaries, matching plink's `--split-x b38`.
X_PAR1 = (10_001, 2_781_479)
X_NONPAR = (2_781_480, 155_701_382)
X_PAR2 = (155_701_383, 156_030_895)
X_PAR1_N = 30
X_NONPAR_N = 250
X_PAR2_N = 20

# Male-specific region of chrY (GRCh38), away from PAR1/PAR2.
Y_REGION = (2_786_855, 26_600_000)
Y_N = 30

# Extra dense block placed inside the chr6 MHC high-LD interval.
MHC_BLOCK = (25_400_000, 33_400_000)
MHC_BLOCK_N = 30

FST_SUPER = 0.12
FST_SUB = 0.01
ANC_FREQ_RANGE = (0.10, 0.50)
FREQ_CLAMP = (0.03, 0.97)

SUPER_POPS = ["AFR", "AMR", "CSA", "EAS", "EUR"]
SUB_POPS = {
    "AFR": ["YRI", "ESN"],
    "AMR": ["MXL", "PEL"],
    "CSA": ["GIH", "PJL"],
    "EAS": ["CHB", "JPT"],
    "EUR": ["GBR", "TSI"],
}
REF_PER_SUBPOP = 25  # -> 50 reference samples per super-population, 250 total

STRATUM_N = 125
STRATA = [
    {"id": "study_eur", "ancestry": "EUR", "prefix": "EUR"},
    {"id": "study_afr", "ancestry": "AFR", "prefix": "AFR"},
]

# Defect slots (1-based sample number inside a stratum).  Sample IDs are
# "<PREFIX>_S<slot:03d>", so the slot fully determines the ID.
SLOT_SEX_FAMFEMALE_GENOMALE = 10
SLOT_SEX_FAMMALE_GENOFEMALE = 11
SLOT_MZ_A, SLOT_MZ_B = 20, 21
SLOT_PARENT, SLOT_OFFSPRING = 30, 31
SLOT_SIB_A, SLOT_SIB_B = 32, 33
SLOT_MISS_SAMPLE = {40: 0.15, 41: 0.22, 42: 0.05}
SLOT_ANCESTRY_OUTLIER = {
    "study_eur": {50: "EAS", 51: "AFR", 52: "AMR"},
    "study_afr": {50: "EUR", 51: "CSA", 52: "EAS"},
}
# Deterministic carriers for the planted rare variants.
SLOT_RARE_CARRIERS = (60, 61)

# Variant defect blocks: (chrom, minimum position, how many).
BLOCK_VARMISS_HIGH = (3, 120_000_000, 20)
BLOCK_VARMISS_BORDER = (4, 100_000_000, 5)
BLOCK_HWE = (5, 150_000_000, 5)
BLOCK_RARE = (6, 100_000_000, 6)
VARMISS_HIGH_RANGE = (0.12, 0.25)
VARMISS_BORDER_RANGE = (0.04, 0.06)

FLIP_FRACTION = 0.015      # fraction of autosomal variants strand-flipped
N_AMBIGUOUS = 40           # autosomal variants given A/T or C/G alleles

# GRCh38 long-range / high-LD regions.  Identical interval set to the legacy
# GLAD/EDGI/NBR asset `highLDregions4bim_b38.awk` (1-based inclusive there),
# re-expressed here as 0-based half-open BED.
HIGH_LD_REGIONS_1BASED = [
    (1, 47822309, 51822307, "hld_1_1"),
    (2, 85861220, 100425020, "hld_2_1"),
    (2, 133908698, 137408698, "hld_2_2_LCT"),
    (2, 182309768, 189309768, "hld_2_3"),
    (3, 47483507, 49987563, "hld_3_1"),
    (3, 83368160, 86868160, "hld_3_2"),
    (3, 88868161, 96298466, "hld_3_3"),
    (5, 44464142, 51168409, "hld_5_1"),
    (5, 98636397, 101136397, "hld_5_2"),
    (5, 129636409, 132636409, "hld_5_3"),
    (5, 136136413, 139136412, "hld_5_4"),
    (6, 25391794, 33424245, "hld_6_1_MHC"),
    (6, 57027244, 63232136, "hld_6_2"),
    (6, 139637171, 142137170, "hld_6_3"),
    (7, 55158099, 67090863, "hld_7_1"),
    (8, 8105069, 12105082, "hld_8_1_INV8p23"),
    (8, 43025701, 48924888, "hld_8_2"),
    (8, 110918596, 113918595, "hld_8_3"),
    (10, 36671067, 43184546, "hld_10_1"),
    (11, 46021874, 57475951, "hld_11_1"),
    (11, 88127185, 91127184, "hld_11_2"),
    (12, 32955800, 41319931, "hld_12_1"),
    (12, 110599476, 113099475, "hld_12_2"),
    (20, 33948534, 36438183, "hld_20_1"),
]

UNAMBIGUOUS_PAIRS = [
    ("A", "G"), ("G", "A"), ("C", "T"), ("T", "C"),
    ("A", "C"), ("C", "A"), ("G", "T"), ("T", "G"),
]
AMBIGUOUS_PAIRS = [("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")]
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

MISSING = 3
# PLINK 1 .bed two-bit codes.  00 = homozygous for the first .bim allele (A1),
# 01 = missing, 10 = heterozygous, 11 = homozygous for the second allele (A2).
# Genotypes are stored here as the A1 dosage.
BED_CODE = {2: 0b00, 1: 0b10, 0: 0b11, MISSING: 0b01}


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #

def clamp(value, low, high):
    return low if value < low else (high if value > high else value)


def bn_draw(rng, p, fst):
    """Balding-Nichols: draw a population frequency around ancestral ``p``."""
    if fst <= 0:
        return p
    scale = (1.0 - fst) / fst
    return rng.betavariate(p * scale, (1.0 - p) * scale)


def draw_diploid(rng, p):
    return (1 if rng.random() < p else 0) + (1 if rng.random() < p else 0)


def draw_haploid_hom(rng, p):
    """Hemizygous call, written as a homozygote (PLINK's haploid convention)."""
    return 2 if rng.random() < p else 0


def transmit(rng, parent_dosage, p):
    """One allele transmitted from ``parent_dosage``, one from the population."""
    from_parent = 1 if rng.random() < parent_dosage / 2.0 else 0
    from_pool = 1 if rng.random() < p else 0
    return from_parent + from_pool


def spread_positions(rng, low, high, count):
    """``count`` distinct sorted positions in [low, high]."""
    span = high - low + 1
    if count > span:
        raise ValueError("more variants requested than positions available")
    picked = set()
    while len(picked) < count:
        picked.add(rng.randint(low, high))
    return sorted(picked)


CHROM_ORDER = {str(c): c for c in range(1, 23)}
CHROM_ORDER["X"] = 23
CHROM_ORDER["Y"] = 24


# --------------------------------------------------------------------------- #
# Variant panel
# --------------------------------------------------------------------------- #

class Variant:
    __slots__ = ("chrom", "vid", "pos", "a1", "a2", "region", "freq", "ambiguous", "flipped")

    def __init__(self, chrom, vid, pos, region):
        self.chrom = chrom
        self.vid = vid
        self.pos = pos
        self.region = region      # "auto" | "par1" | "nonpar" | "par2" | "Y"
        self.a1 = "A"
        self.a2 = "G"
        self.freq = {}            # super-pop -> A1 frequency
        self.ambiguous = False
        self.flipped = False


def build_variants(seed):
    rng = random.Random(seed + 101)

    total_len = sum(CHROM_LEN_B38.values())
    per_chrom = {}
    for chrom, length in CHROM_LEN_B38.items():
        n = int(round(AUTOSOME_TARGET * length / total_len))
        per_chrom[chrom] = max(AUTOSOME_MIN_PER_CHROM, n)

    variants = []
    for chrom in range(1, 23):
        low = VARIANT_MARGIN
        high = CHROM_LEN_B38[chrom] - VARIANT_MARGIN
        positions = set(spread_positions(rng, low, high, per_chrom[chrom]))
        if chrom == 6:
            # Guarantee a dense block inside the chr6 MHC high-LD interval.
            positions |= set(spread_positions(rng, MHC_BLOCK[0], MHC_BLOCK[1], MHC_BLOCK_N))
        for pos in sorted(positions):
            variants.append(Variant(str(chrom), None, pos, "auto"))

    for low, high, count, region in (
        (X_PAR1[0] + 10_000, X_PAR1[1], X_PAR1_N, "par1"),
        (X_NONPAR[0], X_NONPAR[1], X_NONPAR_N, "nonpar"),
        (X_PAR2[0], X_PAR2[1], X_PAR2_N, "par2"),
    ):
        for pos in spread_positions(rng, low, high, count):
            variants.append(Variant("X", None, pos, region))

    for pos in spread_positions(rng, Y_REGION[0], Y_REGION[1], Y_N):
        variants.append(Variant("Y", None, pos, "Y"))

    variants.sort(key=lambda v: (CHROM_ORDER[v.chrom], v.pos))

    # rsID-style identifiers: strictly increasing, pseudo-random increments.
    rs = 1_000_003
    for variant in variants:
        rs += rng.randint(3, 997)
        variant.vid = "rs%d" % rs

    # Alleles and Balding-Nichols frequencies.
    freq_rng = random.Random(seed + 202)
    for variant in variants:
        variant.a1, variant.a2 = UNAMBIGUOUS_PAIRS[freq_rng.randrange(len(UNAMBIGUOUS_PAIRS))]
        anc = freq_rng.uniform(*ANC_FREQ_RANGE)
        for pop in SUPER_POPS:
            variant.freq[pop] = clamp(bn_draw(freq_rng, anc, FST_SUPER), *FREQ_CLAMP)
        for pop in SUPER_POPS:
            for sub in SUB_POPS[pop]:
                variant.freq[sub] = clamp(
                    bn_draw(freq_rng, variant.freq[pop], FST_SUB), *FREQ_CLAMP
                )
    return variants


def pick_run(variants, chrom, min_pos, count):
    """Indices of the first ``count`` variants on ``chrom`` at/after ``min_pos``."""
    out = [
        i for i, v in enumerate(variants)
        if v.chrom == str(chrom) and v.pos >= min_pos
    ]
    if len(out) < count:
        raise ValueError("not enough variants on chr%s after %d" % (chrom, min_pos))
    return out[:count]


def in_high_ld(chrom, pos):
    for hc, start, end, _name in HIGH_LD_REGIONS_1BASED:
        if str(hc) == chrom and start <= pos <= end:
            return True
    return False


# --------------------------------------------------------------------------- #
# Samples
# --------------------------------------------------------------------------- #

class Sample:
    __slots__ = ("fid", "iid", "fam_sex", "bio_sex", "pop", "note")

    def __init__(self, iid, fam_sex, bio_sex, pop, note=""):
        self.fid = iid
        self.iid = iid
        self.fam_sex = fam_sex    # what the .fam records: 1 male, 2 female
        self.bio_sex = bio_sex    # what the genotypes look like
        self.pop = pop
        self.note = note


def build_reference_samples():
    samples = []
    meta = []
    for super_pop in SUPER_POPS:
        for sub in SUB_POPS[super_pop]:
            for i in range(1, REF_PER_SUBPOP + 1):
                iid = "%s_%03d" % (sub, i)
                sex = 1 if i % 2 else 2
                samples.append(Sample(iid, sex, sex, sub))
                meta.append((iid, iid, sub, super_pop))
    return samples, meta


def build_study_samples(stratum):
    prefix = stratum["prefix"]
    home_pop = stratum["ancestry"]
    outliers = SLOT_ANCESTRY_OUTLIER[stratum["id"]]
    samples = []
    for slot in range(1, STRATUM_N + 1):
        iid = "%s_S%03d" % (prefix, slot)
        sex = 1 if slot % 2 else 2
        fam_sex, bio_sex = sex, sex
        pop = home_pop
        note = ""

        if slot == SLOT_SEX_FAMFEMALE_GENOMALE:
            fam_sex, bio_sex, note = 2, 1, "sex_mismatch_fam_female_geno_male"
        elif slot == SLOT_SEX_FAMMALE_GENOFEMALE:
            fam_sex, bio_sex, note = 1, 2, "sex_mismatch_fam_male_geno_female"
        elif slot == SLOT_MZ_A:
            fam_sex = bio_sex = 2
            note = "mz_pair_donor"
        elif slot == SLOT_MZ_B:
            fam_sex = bio_sex = 2
            note = "mz_pair_copy"
        elif slot == SLOT_PARENT:
            fam_sex = bio_sex = 1
            note = "parent_of_%s_S%03d" % (prefix, SLOT_OFFSPRING)
        elif slot == SLOT_OFFSPRING:
            fam_sex = bio_sex = 2
            note = "offspring_of_%s_S%03d" % (prefix, SLOT_PARENT)
        elif slot in (SLOT_SIB_A, SLOT_SIB_B):
            fam_sex = bio_sex = 1 if slot == SLOT_SIB_A else 2
            note = "full_sib"
        elif slot in SLOT_MISS_SAMPLE:
            # Forced male so that chrY contributes no baseline missingness and
            # the achieved F_MISS equals the planted target exactly.
            fam_sex = bio_sex = 1
            note = "sample_missingness_%.2f" % SLOT_MISS_SAMPLE[slot]
        elif slot in outliers:
            pop = outliers[slot]
            note = "ancestry_outlier_%s" % pop

        samples.append(Sample(iid, fam_sex, bio_sex, pop, note))
    return samples


# --------------------------------------------------------------------------- #
# Genotype generation
# --------------------------------------------------------------------------- #

def draw_column(rng, variant, samples, pop_override=None):
    """One variant's genotypes for independent samples (derived ones skipped)."""
    col = bytearray(len(samples))
    for idx, sample in enumerate(samples):
        pop = pop_override or sample.pop
        p = variant.freq[pop]
        if variant.region == "Y":
            col[idx] = draw_haploid_hom(rng, p) if sample.bio_sex == 1 else MISSING
        elif variant.region == "nonpar" and sample.bio_sex == 1:
            col[idx] = draw_haploid_hom(rng, p)
        else:
            col[idx] = draw_diploid(rng, p)
    return col


def generate_reference(variants, samples, seed):
    rng = random.Random(seed + 303)
    autosomal = [v for v in variants if v.region == "auto"]
    geno = []
    for variant in autosomal:
        col = bytearray(len(samples))
        for idx, sample in enumerate(samples):
            col[idx] = draw_diploid(rng, variant.freq[sample.pop])
        geno.append(col)
    return autosomal, geno


def generate_stratum(variants, samples, stratum, seed, defects):
    rng = random.Random(seed + 404 + sum(ord(c) for c in stratum["id"]))
    home_pop = stratum["ancestry"]
    index = {s.iid: i for i, s in enumerate(samples)}
    prefix = stratum["prefix"]

    def slot_index(slot):
        return index["%s_S%03d" % (prefix, slot)]

    i_mz_a, i_mz_b = slot_index(SLOT_MZ_A), slot_index(SLOT_MZ_B)
    i_par, i_off = slot_index(SLOT_PARENT), slot_index(SLOT_OFFSPRING)
    i_sib_a, i_sib_b = slot_index(SLOT_SIB_A), slot_index(SLOT_SIB_B)
    derived = {i_mz_b, i_off, i_sib_a, i_sib_b}

    geno = []
    for variant in variants:
        col = bytearray(len(samples))
        for idx, sample in enumerate(samples):
            if idx in derived:
                continue
            p = variant.freq[sample.pop]
            if variant.region == "Y":
                col[idx] = draw_haploid_hom(rng, p) if sample.bio_sex == 1 else MISSING
            elif variant.region == "nonpar" and sample.bio_sex == 1:
                col[idx] = draw_haploid_hom(rng, p)
            else:
                col[idx] = draw_diploid(rng, p)

        # Monozygotic / duplicate pair: an exact genotype copy.
        col[i_mz_b] = col[i_mz_a]

        p_home = variant.freq[home_pop]
        if variant.region == "auto" or variant.region in ("par1", "par2"):
            # Parent-offspring: one allele transmitted, one drawn from the pool.
            col[i_off] = transmit(rng, col[i_par], p_home)
            # Full sibs: two virtual (not-genotyped) parents, both sibs
            # inherit one allele from each.
            vp1 = draw_diploid(rng, p_home)
            vp2 = draw_diploid(rng, p_home)
            for i_sib in (i_sib_a, i_sib_b):
                a1 = 1 if rng.random() < vp1 / 2.0 else 0
                a2 = 1 if rng.random() < vp2 / 2.0 else 0
                col[i_sib] = a1 + a2
        else:
            # chrX non-PAR and chrY: drawn independently for the derived
            # samples too (KING relatedness uses autosomes).
            for i_sib in (i_off, i_sib_a, i_sib_b):
                sample = samples[i_sib]
                if variant.region == "Y":
                    col[i_sib] = draw_haploid_hom(rng, p_home) if sample.bio_sex == 1 else MISSING
                elif sample.bio_sex == 1:
                    col[i_sib] = draw_haploid_hom(rng, p_home)
                else:
                    col[i_sib] = draw_diploid(rng, p_home)
        geno.append(col)

    # ---- planted variant-level defects -----------------------------------
    for vi in defects["hwe"]:
        for idx in range(len(samples)):
            geno[vi][idx] = 1  # every sample heterozygous

    for order, vi in enumerate(defects["rare"]):
        n_carriers = 1 if order < len(defects["rare"]) // 2 else 2
        for idx in range(len(samples)):
            geno[vi][idx] = 0
        for slot in SLOT_RARE_CARRIERS[:n_carriers]:
            geno[vi][slot_index(slot)] = 1

    miss_rng = random.Random(seed + 505 + sum(ord(c) for c in stratum["id"]))
    var_miss_rates = {}
    for vi in defects["varmiss_high"]:
        var_miss_rates[vi] = miss_rng.uniform(*VARMISS_HIGH_RANGE)
    for vi in defects["varmiss_border"]:
        var_miss_rates[vi] = miss_rng.uniform(*VARMISS_BORDER_RANGE)
    for vi, rate in sorted(var_miss_rates.items()):
        n_missing = int(round(rate * len(samples)))
        for idx in miss_rng.sample(range(len(samples)), n_missing):
            geno[vi][idx] = MISSING

    # ---- planted sample-level missingness --------------------------------
    n_var = len(variants)
    for slot, target in sorted(SLOT_MISS_SAMPLE.items()):
        idx = slot_index(slot)
        target_n = int(round(target * n_var))
        callable_vi = [vi for vi in range(n_var) if geno[vi][idx] != MISSING]
        already = n_var - len(callable_vi)
        extra = target_n - already
        if extra > 0:
            for vi in miss_rng.sample(callable_vi, extra):
                geno[vi][idx] = MISSING
    return geno


# --------------------------------------------------------------------------- #
# Writers
# --------------------------------------------------------------------------- #

def write_bed(path, geno, n_samples):
    n_bytes = (n_samples + 3) // 4
    with open(path, "wb") as handle:
        handle.write(bytes([0x6C, 0x1B, 0x01]))
        for col in geno:
            buf = bytearray(n_bytes)
            for idx, dosage in enumerate(col):
                buf[idx >> 2] |= BED_CODE[dosage] << ((idx & 3) * 2)
            handle.write(bytes(buf))


def write_bim(path, variants, flip=False):
    with open(path, "w", newline="\n") as handle:
        for v in variants:
            a1, a2 = v.a1, v.a2
            if flip and v.flipped:
                a1, a2 = COMPLEMENT[a1], COMPLEMENT[a2]
            handle.write("%s\t%s\t0\t%d\t%s\t%s\n" % (v.chrom, v.vid, v.pos, a1, a2))


def write_fam(path, samples):
    with open(path, "w", newline="\n") as handle:
        for s in samples:
            handle.write("%s %s 0 0 %d -9\n" % (s.fid, s.iid, s.fam_sex))


def write_lines(path, lines):
    with open(path, "w", newline="\n") as handle:
        for line in lines:
            handle.write(line + "\n")


# --------------------------------------------------------------------------- #
# Statistics (for the README)
# --------------------------------------------------------------------------- #

def variant_stats(geno_col):
    called = [g for g in geno_col if g != MISSING]
    n = len(called)
    f_miss = (len(geno_col) - n) / len(geno_col)
    if n == 0:
        return 0.0, f_miss, 0
    a1 = sum(called)
    freq = a1 / (2.0 * n)
    maf = min(freq, 1.0 - freq)
    n_het = sum(1 for g in called if g == 1)
    return maf, f_miss, n_het


def sample_f_miss(geno, idx):
    missing = sum(1 for col in geno if col[idx] == MISSING)
    return missing / len(geno)


# --------------------------------------------------------------------------- #
# README
# --------------------------------------------------------------------------- #

# Recorded output of the verification runs described below.  Kept in the
# generator so that a plain `python3 scripts/gwasqc/generate_fixtures.py`
# reproduces data/gwasqc/README.md byte-for-byte.  `--verification FILE`
# replaces it when the checks are re-run.
RECORDED_VERIFICATION = r"""
Every check below passed.  Commands are written as if `plink2` and `plink` were
on PATH; they were run inside the two containers above, from this directory, and
re-running them on the committed data reproduces the numbers quoted.

### 1. All three filesets load, and `--freq` / `--missing` / `--hardy` run

```
plink2 --bfile genotypes/study_eur   --freq --missing --hardy --out out_eur
plink2 --bfile genotypes/study_afr   --freq --missing --hardy --out out_afr
plink2 --bfile reference/ref_panel   --freq --missing --hardy --out out_ref
```

| Fileset | Loaded |
| --- | --- |
| `genotypes/study_eur` | `125 samples (61 females, 64 males)`, `3962 variants` |
| `genotypes/study_afr` | `125 samples (61 females, 64 males)`, `3962 variants` |
| `reference/ref_panel` | `250 samples (120 females, 130 males)`, `3632 variants` |

No errors. `--hardy` reports `Skipping 30 haploid variants` (chrY) and writes a
separate `.hardy.x` for chrX, as expected. `plink2 --bfile ... --make-pgen`
followed by `plink2 --pfile ...` round-trips cleanly, so the pgen route works
too.

### 2. Sex check -- both defect a1 and a2 are caught, and nothing else

plink1.9 route (`--split-x` + `--check-sex ycount`):

```
plink  --bfile genotypes/study_eur --split-x b38 no-fail --make-bed --out sx_eur
plink  --bfile sx_eur --check-sex ycount --out cs_eur
```

`--split-x: 50 chromosome codes changed` (exactly the 30 PAR1 + 20 PAR2
variants), then
`--check-sex: 250 Xchr and 30 Ychr variant(s) scanned, 2 problems detected`:

```
FID       IID       PEDSEX SNPSEX STATUS  F        YCOUNT
EUR_S010  EUR_S010  2      1      PROBLEM 1        30
EUR_S011  EUR_S011  1      2      PROBLEM 0.01474  0
AFR_S010  AFR_S010  2      1      PROBLEM 1        30
AFR_S011  AFR_S011  1      2      PROBLEM 0.1221   0
```

F separates cleanly by genotype sex: `PEDSEX 1` mean F 0.985 / 0.986, `PEDSEX 2`
mean F 0.011 / -0.005; YCOUNT is 29.8-29.9 for `SNPSEX 1` and exactly 0.0 for
`SNPSEX 2`.

plink2 route -- **plink2 2.0.0a.6.9 can do the whole sex check**, so plink1.9 is
not required for this stage.  The flags differ: it is `--split-par`, not
`--split-x`, and `--check-sex` takes explicit thresholds:

```
plink2 --bfile genotypes/study_eur --split-par b38 --make-bed --out p2sx_eur
plink2 --bfile p2sx_eur --check-sex max-female-xf=0.2 min-male-xf=0.8 \
       max-female-ycount=0 min-male-ycount=1 --out p2cs_eur
```

`--split-par: 50 chromosome codes changed`, then
`--check-sex: 250 chrX variants and 30 variants scanned, 2 problems detected` --
the same two samples per stratum (`F` 1 / 0.0158 for EUR, 1 / 0.1228 for AFR;
`YRATE` 1 / 0).  `plink2 --split-x` is an `Unrecognized flag`, and plink2
`--check-sex` with no thresholds at all uses `min-male-xf=1`,
`max-female-yrate=0`, which on unsplit data flags 65 of 125 samples -- so both
the PAR split and explicit thresholds matter.

### 3. Relatedness -- KING finds all three planted pairs and no false positives

```
plink2 --bfile genotypes/study_eur --chr 1-22 \
       --make-king-table --king-table-filter 0.08 --out king_eur
```

Complete output (every pair above kinship 0.08):

```
IID1      IID2      NSNP  HETHET    IBS0       KINSHIP
EUR_S021  EUR_S020  3620  0.348066  0          0.5
EUR_S031  EUR_S030  3627  0.175076  0          0.251235
EUR_S033  EUR_S032  3624  0.204194  0.0140728  0.258607
AFR_S021  AFR_S020  3629  0.336181  0          0.5
AFR_S031  AFR_S030  3628  0.168964  0          0.243755
AFR_S033  AFR_S032  3625  0.216552  0.0184828  0.252937
```

Kinship 0.5 with IBS0 = 0 for the duplicate/MZ pair; ~0.25 with IBS0 = 0 for
parent-offspring; ~0.25 with IBS0 > 0 for the full sibs -- i.e. the two
first-degree classes are distinguishable by IBS0, as KING's decision rule
requires.  Nothing else in either stratum exceeds 0.08.

### 4. Missingness

From `out_*.smiss` / `out_*.vmiss`.  The three planted samples are the top three
by F_MISS in both strata (0.2201, 0.1499, 0.0500); the next highest sample is
0.0083.  Exactly 20 variants exceed F_MISS 0.10 in each stratum, all on chr3 and
all from block d1 (range 0.128-0.256); 11 variants sit in (0.03, 0.10], the 5
planted chr4 borderline variants (0.048-0.064) plus 6 that get there by chance
from the sample-level missingness.

### 5. HWE

Exactly 5 autosomal variants have `--hardy` p < 1e-10 in each stratum, and they
are precisely the 5 planted chr5 variants (p = 4.7e-37 to 1.9e-36, `O(HET_A1)`
= 1, `E(HET_A1)` = 0.5).  The next most extreme unplanted variant has
p = 5.2e-05, so any threshold between those two isolates the plants.

### 6. MAF

The six lowest `ALT_FREQS` values in each stratum are the six planted chr6
variants (0.0040 / 0.0040 / 0.0040 and 0.0081 / 0.0081 / 0.0081 in `study_eur`).
Counting autosomal variants with MAF < 0.01: 8 in `study_eur`, 7 in
`study_afr`, 0 in the reference panel -- the 6 plants plus 2 (resp. 1) that fall
below by sampling chance.  Reference-panel mean MAF is 0.292.

### 7. High-LD overlap

```
awk 'NR==FNR {n++; c[n]=$1; s[n]=$2; e[n]=$3; nm[n]=$4; next}
     {for (i=1;i<=n;i++) if ($1==c[i] && $4>s[i] && $4<=e[i]) {hit++; cnt[nm[i]]++; next}}
     END {print hit; for (k in cnt) print k, cnt[k]}' \
  highld/high_ld_regions_b38.bed genotypes/study_eur.bim
```

202 of the 3632 autosomal variants fall inside a listed interval, spread over
**all 24** intervals (2-38 variants each); 38 of them are in `hld_6_1_MHC`.  The
`$4 > start && $4 <= end` comparison is the correct one for a 0-based half-open
BED against a 1-based `.bim` position.

### 8. Harmonisation -- the strand flips are exactly the documented set

A naive merge without harmonisation fails, and fails on precisely the right
variants:

```
plink2 --bfile genotypes/study_eur --chr 1-22 --make-bed --out eur_auto
plink  --bfile reference/ref_panel --bmerge eur_auto --make-bed --out merged_naive
# -> Error: 54 variants with 3+ alleles present.  merged_naive-merge.missnp
```

`merged_naive-merge.missnp` has 54 lines and is **set-identical** to the flipped
variant list derived independently by diffing the reference and study `.bim`
allele columns (and to the ID list in the appendix below).  Flipping them makes
the merge succeed:

```
plink  --bfile eur_auto --flip merged_naive-merge.missnp --make-bed --out eur_flip
plink  --bfile afr_auto --flip merged_naive-merge.missnp --make-bed --out afr_flip
plink  --bfile reference/ref_panel --bmerge eur_flip --make-bed --out m1   # 375 people
plink  --bfile m1 --bmerge afr_flip --make-bed --out merged                # 500 people, 3632 variants
```

### 9. Ancestry -- PCA separates the super-populations and isolates the outliers

```
plink2 --bfile merged --maf 0.05 --pca 6 --out pca
```

Eigenvalues 40.4, 22.5, 14.2, 13.9, 2.1, 2.0 -- four structure PCs for five
populations, then a drop of ~7x, which is what five discrete populations should
give.  Nearest-reference-centroid assignment on PC1-PC4:

* reference self-assignment **250/250 correct**;
* minimum between-centroid distance 0.107 against a mean within-cluster
  distance of 0.007 (max 0.016), i.e. clusters separated by ~15x their own
  radius;
* `study_eur`: 122/125 assigned EUR; the 3 exceptions are exactly `EUR_S050`
  -> EAS, `EUR_S051` -> AFR, `EUR_S052` -> AMR;
* `study_afr`: 122/125 assigned AFR; the 3 exceptions are exactly `AFR_S050`
  -> EUR, `AFR_S051` -> CSA, `AFR_S052` -> EAS.

Every planted ancestry outlier lands on the super-population it was drawn
from, and no unplanted sample is misassigned.

### 10. Remaining pipeline stages smoke-tested

```
plink2 --bfile genotypes/study_eur --chr 1-22 --het --out het_eur
plink2 --bfile genotypes/study_eur --chr 1-22 --indep-pairwise 1500 150 0.2 --out ld_eur
plink2 --bfile reference/ref_panel --keep reference/keep/EUR.keep --freq --out k_EUR
```

`--het` runs (mean F = -0.0009 over 125 samples, as it should be for
HWE-conforming data).  `--indep-pairwise 1500 150 0.2` removes only 15 of 3632
variants -- see the caveat about there being no LD structure.  Each of the five
`keep/<SUPERPOP>.keep` files yields `--keep: 50 samples remaining`.

### Not verified

* **No Nextflow / nf-test run.** These fixtures have not yet been consumed by
  the pipeline; only the PLINK-level behaviour above is checked.
* **Heterozygosity outliers are not planted**, so there is nothing for the het
  stage to *find*; only that `--het` runs and produces sane values.
* **Ancestry was verified by merge + PCA, not by the pipeline's
  `--pca allele-wts` + `--score` projection route.** Separability is
  established; the projection path itself is untested here.
* **LD pruning is not meaningfully exercised** (no LD structure by
  construction), and consequently neither is the interaction between pruning
  and the high-LD exclusion beyond coordinate overlap.
* **No MAF/HWE/missingness defects in the reference panel**, so a pipeline that
  QCs the reference has nothing to catch there.
"""

VERIFICATION = """\
## Verification

Verified with the pinned containers

* `quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0` (PLINK v2.0.0-a.6.9LM, 29 Jan 2025)
* `quay.io/biocontainers/plink:1.90b6.21--h779adbc_1` (PLINK v1.90b6.21, 19 Oct 2020)

@@VERIFICATION_BODY@@
"""


def build_readme(ctx):
    L = []
    A = L.append
    A("# gwasqc test fixtures")
    A("")
    A("Synthetic genotype fixtures for the `gwasqc` genotype-QC pipeline "
      "(<https://github.com/lyh970817/gwasqc>).")
    A("**No real human data.** Every genotype here is simulated.")
    A("")
    A("Everything in this directory is reproduced from scratch by")
    A("")
    A("```bash")
    A("python3 scripts/gwasqc/generate_fixtures.py   # seed %d, stdlib only" % ctx["seed"])
    A("```")
    A("")
    A("The generator is deterministic: a re-run with the same seed reproduces "
      "byte-identical files. Changing the script changes the bytes, so treat "
      "the committed data as the artefact and the script as its provenance.")
    A("")

    # -------------------------------------------------- inventory
    A("## Files")
    A("")
    A("| Path | Size | Contents |")
    A("| --- | --- | --- |")
    for rel, desc in ctx["inventory"]:
        size = ctx["sizes"][rel]
        A("| `%s` | %s | %s |" % (rel, human(size), desc))
    A("")
    A("Total, excluding this README: **%s**." % human(sum(ctx["sizes"].values())))
    A("")

    # -------------------------------------------------- conventions
    A("## Conventions")
    A("")
    A("* **Build: GRCh38 / b38.** Positions are GRCh38 primary-assembly "
      "coordinates; chrX PAR boundaries match plink1.9 `--split-x b38` and "
      "plink2 `--split-par b38` "
      "(PAR1 `X:%d-%d`, PAR2 `X:%d-%d`)."
      % (X_PAR1[0], X_PAR1[1], X_PAR2[0], X_PAR2[1]))
    A("* **Chromosome codes** in `.bim` are bare PLINK codes without a `chr` "
      "prefix: `1`-`22`, `X`, `Y` (plus `XY` after a PAR split). "
      "`high_ld_regions_b38.bed` uses the same "
      "codes so it can be joined straight onto a `.bim`/`.pvar`.")
    A("* `.bim` is tab-delimited, `.fam` space-delimited (what PLINK itself "
      "writes). `.bim` column 5 is A1 (the counted allele), column 6 is A2.")
    A("* `.fam`/`ref_panel.fam` phenotype column is `-9` throughout: this is a "
      "QC fixture, not an association fixture.")
    A("* `FID == IID` for every sample, and the `.fam` PAT/MAT columns are `0`. "
      "The related pairs below carry **no pedigree**: relatedness has to be "
      "inferred from genotypes, which is the point.")
    A("* `high_ld_regions_b38.bed` is a real BED: **0-based, half-open**. Its "
      "interval set is identical to the legacy asset "
      "`GLAD_EDGI_NBR_qc_workflows/full_pipe/dependencies/highLDregions4bim_b38.awk` "
      "(which states the same intervals 1-based inclusive).")
    A("")

    # -------------------------------------------------- composition
    A("## Composition")
    A("")
    A("### Variants (%d in the study filesets, %d in the reference)"
      % (ctx["n_var"], ctx["n_var_auto"]))
    A("")
    A("| Region | Variants |")
    A("| --- | --- |")
    A("| autosomes chr1-22 | %d |" % ctx["n_var_auto"])
    A("| chrX PAR1 (`X:%d-%d`) | %d |" % (X_PAR1[0], X_PAR1[1], ctx["n_par1"]))
    A("| chrX non-PAR | %d |" % ctx["n_nonpar"])
    A("| chrX PAR2 (`X:%d-%d`) | %d |" % (X_PAR2[0], X_PAR2[1], ctx["n_par2"]))
    A("| chrY (male-specific region) | %d |" % ctx["n_y"])
    A("")
    A("Per-autosome counts: " + ", ".join(
        "chr%s %d" % (c, n) for c, n in ctx["per_chrom"]) + ".")
    A("")
    A("Variant IDs are rsID-style and strictly increasing with genomic order. "
      "The %d chrX PAR variants exist so that the PAR split has work to do: "
      "both plink1.9 `--split-x b38` and plink2 `--split-par b38` move exactly "
      "these %d variants to `XY`/`25`, leaving the %d non-PAR chrX variants "
      "that the sex check uses. plink2 has **no** `--split-x`; there the flag "
      "is spelled `--split-par`."
      % (ctx["n_par1"] + ctx["n_par2"], ctx["n_par1"] + ctx["n_par2"],
         ctx["n_nonpar"]))
    A("")
    A("### Samples")
    A("")
    A("| Fileset | Samples | `.fam` sex 1 / 2 |")
    A("| --- | --- | --- |")
    for name, n, nm, nf in ctx["sample_summary"]:
        A("| `%s` | %d | %d / %d |" % (name, n, nm, nf))
    A("")
    A("Study sample IDs are `<PREFIX>_S<nnn>` with `PREFIX` = `EUR`/`AFR` and "
      "`nnn` = 001-%d; reference sample IDs are `<SUBPOP>_<nnn>`. The numeric "
      "part of a study ID is stable across regenerations, so the defect IDs "
      "below never move." % STRATUM_N)
    A("")

    # -------------------------------------------------- samplesheet
    A("### `samplesheet.csv`")
    A("")
    A("```csv")
    for line in ctx["samplesheet"]:
        A(line)
    A("```")
    A("")
    A("Header is exactly `id,ancestry,bed,bim,fam`; the three path columns are "
      "**relative to this directory** (`data/gwasqc/`).")
    A("")

    # -------------------------------------------------- ancestry
    A("## Ancestry design")
    A("")
    A("Allele frequencies come from a two-level Balding-Nichols model: one "
      "ancestral frequency per variant drawn from U(%.2f, %.2f), differentiated "
      "across the five super-populations with **Fst = %.2f**, then across two "
      "sub-populations per super-population with Fst = %.2f. Frequencies are "
      "clamped to [%.2f, %.2f]. Genotypes are binomial draws, so every variant "
      "is HWE-conforming unless a defect was planted into it."
      % (ANC_FREQ_RANGE[0], ANC_FREQ_RANGE[1], FST_SUPER, FST_SUB,
         FREQ_CLAMP[0], FREQ_CLAMP[1]))
    A("")
    A("| Super-population | Sub-populations | Reference samples |")
    A("| --- | --- | --- |")
    for pop in SUPER_POPS:
        A("| %s | %s | %d |" % (pop, ", ".join(SUB_POPS[pop]), 2 * REF_PER_SUBPOP))
    A("")
    A("`reference/ref_pop.tsv` is tab-separated with the header "
      "`FID IID pop super_pop`; `reference/keep/<SUPERPOP>.keep` lists "
      "`FID IID` (space-separated, no header) for that super-population.")
    A("")
    A("Study strata are drawn from their matching super-population frequencies "
      "(`study_eur` from EUR, `study_afr` from AFR), so a reference-projected "
      "PCA should place them on top of the matching reference cluster -- except "
      "for the planted ancestry outliers.")
    A("")
    A("**Reference/study harmonisation.** Reference and study filesets share "
      "variant ID, position and alleles on all %d autosomal variants, with two "
      "documented exceptions that exist to exercise harmonisation:" % ctx["n_var_auto"])
    A("")
    A("* **%d strand-flipped variants** (%.1f%% of autosomal variants): the "
      "study `.bim` carries the complement of the reference alleles "
      "(e.g. reference `A G` -> study `T C`) while the genotype dosages are "
      "unchanged, which is exactly what a strand flip looks like. All of them "
      "have unambiguous (non-`A/T`, non-`C/G`) allele pairs, so the flip is "
      "resolvable. Both strata flip the same variant set."
      % (ctx["n_flipped"], 100.0 * ctx["n_flipped"] / ctx["n_var_auto"]))
    A("* **%d strand-ambiguous variants** carry `A/T` or `C/G` alleles in both "
      "the reference and the study sets (never flipped). These are the ones an "
      "IUPAC/ambiguity check must refuse to resolve." % ctx["n_ambiguous"])
    A("")
    A("Full ID lists: see `flipped variants` and `ambiguous variants` at the "
      "bottom of this file.")
    A("")

    # -------------------------------------------------- defects
    A("## Planted defects")
    A("")
    A("Every ID below is deterministic: the slot number in a study sample ID "
      "encodes its role, and both strata carry the same slots. Sample-level "
      "defect IDs are listed per stratum; variant-level defects hit the same "
      "variant IDs in both strata.")
    A("")
    A("### Sample-level")
    A("")
    A("| # | Defect | IDs | How it was planted | Expected QC outcome |")
    A("| --- | --- | --- | --- | --- |")
    for row in ctx["sample_defects"]:
        A("| " + " | ".join(row) + " |")
    A("")
    A("### Variant-level")
    A("")
    A("| # | Defect | IDs | How it was planted | Expected QC outcome |")
    A("| --- | --- | --- | --- | --- |")
    for row in ctx["variant_defects"]:
        A("| " + " | ".join(row) + " |")
    A("")
    A("### Measured values in the committed data")
    A("")
    A("| Statistic | `study_eur` | `study_afr` |")
    A("| --- | --- | --- |")
    for label, a, b in ctx["measured"]:
        A("| %s | %s | %s |" % (label, a, b))
    A("")

    # -------------------------------------------------- high LD
    A("## High-LD regions")
    A("")
    A("`highld/high_ld_regions_b38.bed` lists %d long-range-LD intervals "
      "(GRCh38), tab-separated, 4 columns `chrom start end name`, 0-based "
      "half-open." % len(HIGH_LD_REGIONS_1BASED))
    A("")
    A("**%d of the %d autosomal variants fall inside a listed interval**, of "
      "which %d sit in the deliberately dense block inside "
      "`hld_6_1_MHC` (`6:%d-%d`). Overlap counts per interval that catches "
      "anything:" % (ctx["n_in_highld"], ctx["n_var_auto"], ctx["n_in_mhc"],
                     MHC_BLOCK[0], MHC_BLOCK[1]))
    A("")
    A("| Interval | Region | Variants |")
    A("| --- | --- | --- |")
    for name, region, count in ctx["highld_overlap"]:
        A("| `%s` | `%s` | %d |" % (name, region, count))
    A("")

    A(ctx["verification"])
    A("")

    # -------------------------------------------------- caveats
    A("## Simplifications and caveats")
    A("")
    A("* The reference panel is **autosomes only** (%d variants). The ancestry "
      "PCA only uses autosomes; nothing here can test a reference chrX path."
      % ctx["n_var_auto"])
    A("* `AMR` and `CSA` are modelled as plain Balding-Nichols populations, "
      "**not** as admixed ones. Real AMR/CSA samples are admixed; a classifier "
      "tuned on this fixture will look better than it is on real data.")
    A("* There is **no LD structure**: variants are drawn independently. A "
      "high-LD region here is high-LD only by coordinate, not by correlation, "
      "so `--indep-pairwise` will prune almost nothing and the high-LD "
      "exclusion asset is exercised by coordinate overlap alone.")
    A("* Defects are planted in the **study strata only**; the reference panel "
      "is clean (no missingness, no HWE violation, no relatives, no sex "
      "mismatch).")
    A("* Relatedness is planted on the **autosomes**; the derived samples' "
      "chrX non-PAR and chrY genotypes are drawn independently.")
    A("* Hemizygous calls (chrX non-PAR and chrY, for samples whose "
      "*genotypes* are male) are written as homozygotes, which is PLINK's "
      "haploid convention. Two PLINK warnings are therefore **expected and "
      "not fixture defects**: before the PAR split every male's diploid call "
      "at the %d PAR variants counts as `het. haploid` (~1.2k calls), and "
      "after the split the ~80 that remain are defect a2's heterozygous chrX "
      "calls; `Nonmissing nonmale Y chromosome genotype(s) present` is "
      "defect a1." % (ctx["n_par1"] + ctx["n_par2"]))
    A("* Planted rare variants are set to a fixed carrier count rather than "
      "sampled, so their MAF is exact. A handful of *unplanted* variants also "
      "fall below MAF 1% by sampling chance; the measured counts above give "
      "the real totals.")
    A("")

    # -------------------------------------------------- ID appendices
    A("## Appendix: variant ID lists")
    A("")
    for title, ids in ctx["id_lists"]:
        A("<details><summary>%s (%d)</summary>" % (title, len(ids)))
        A("")
        A("```")
        for i in range(0, len(ids), 8):
            A(" ".join(ids[i:i + 8]))
        A("```")
        A("")
        A("</details>")
        A("")
    return "\n".join(L) + "\n"


def human(size):
    if size < 1024:
        return "%d B" % size
    if size < 1024 * 1024:
        return "%.1f KB" % (size / 1024.0)
    return "%.2f MB" % (size / (1024.0 * 1024.0))


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="data/gwasqc", type=Path)
    parser.add_argument("--seed", default=DEFAULT_SEED, type=int)
    parser.add_argument("--verification", type=Path, default=None,
                        help="markdown file whose contents become the README's "
                             "verification body (default: the recorded one)")
    args = parser.parse_args()

    out = args.out
    (out / "genotypes").mkdir(parents=True, exist_ok=True)
    (out / "reference" / "keep").mkdir(parents=True, exist_ok=True)
    (out / "highld").mkdir(parents=True, exist_ok=True)

    variants = build_variants(args.seed)
    autosomal_idx = [i for i, v in enumerate(variants) if v.region == "auto"]

    # ---- defect variant selection ---------------------------------------
    defects = {
        "varmiss_high": pick_run(variants, *BLOCK_VARMISS_HIGH),
        "varmiss_border": pick_run(variants, *BLOCK_VARMISS_BORDER),
        "hwe": pick_run(variants, *BLOCK_HWE),
        "rare": pick_run(variants, *BLOCK_RARE),
    }
    reserved = set()
    for vis in defects.values():
        reserved.update(vis)

    pool = [i for i in autosomal_idx if i not in reserved]
    special_rng = random.Random(args.seed + 606)
    n_flip = int(round(FLIP_FRACTION * len(autosomal_idx)))
    chosen = special_rng.sample(pool, n_flip + N_AMBIGUOUS)
    flip_idx = sorted(chosen[:n_flip])
    ambiguous_idx = sorted(chosen[n_flip:])
    for vi in flip_idx:
        variants[vi].flipped = True
    for vi in ambiguous_idx:
        variants[vi].ambiguous = True
        pair = AMBIGUOUS_PAIRS[special_rng.randrange(len(AMBIGUOUS_PAIRS))]
        variants[vi].a1, variants[vi].a2 = pair

    # ---- reference panel -------------------------------------------------
    ref_samples, ref_meta = build_reference_samples()
    ref_variants, ref_geno = generate_reference(variants, ref_samples, args.seed)
    write_bed(out / "reference" / "ref_panel.bed", ref_geno, len(ref_samples))
    write_bim(out / "reference" / "ref_panel.bim", ref_variants, flip=False)
    write_fam(out / "reference" / "ref_panel.fam", ref_samples)
    write_lines(out / "reference" / "ref_pop.tsv",
                ["\t".join(["FID", "IID", "pop", "super_pop"])]
                + ["\t".join(row) for row in ref_meta])
    for super_pop in SUPER_POPS:
        write_lines(out / "reference" / "keep" / ("%s.keep" % super_pop),
                    ["%s %s" % (fid, iid) for fid, iid, _sub, sp in ref_meta
                     if sp == super_pop])

    # ---- study strata ----------------------------------------------------
    stratum_data = {}
    for stratum in STRATA:
        samples = build_study_samples(stratum)
        geno = generate_stratum(variants, samples, stratum, args.seed, defects)
        stem = out / "genotypes" / stratum["id"]
        write_bed(stem.with_suffix(".bed"), geno, len(samples))
        write_bim(stem.with_suffix(".bim"), variants, flip=True)
        write_fam(stem.with_suffix(".fam"), samples)
        stratum_data[stratum["id"]] = (samples, geno)

    # ---- high-LD regions -------------------------------------------------
    write_lines(out / "highld" / "high_ld_regions_b38.bed",
                ["%d\t%d\t%d\t%s" % (chrom, start - 1, end, name)
                 for chrom, start, end, name in HIGH_LD_REGIONS_1BASED])

    # ---- samplesheet -----------------------------------------------------
    samplesheet = ["id,ancestry,bed,bim,fam"] + [
        "%s,%s,genotypes/%s.bed,genotypes/%s.bim,genotypes/%s.fam"
        % (s["id"], s["ancestry"], s["id"], s["id"], s["id"]) for s in STRATA
    ]
    write_lines(out / "samplesheet.csv", samplesheet)

    # ---- README ----------------------------------------------------------
    verification_body = None
    if args.verification and args.verification.exists():
        verification_body = args.verification.read_text()
    readme = build_readme(readme_context(
        args, out, variants, ref_samples, ref_variants, stratum_data,
        defects, flip_idx, ambiguous_idx, samplesheet, verification_body))
    (out / "README.md").write_text(readme)
    print("wrote %s" % (out / "README.md"))


def readme_context(args, out, variants, ref_samples, ref_variants,
                   stratum_data, defects, flip_idx, ambiguous_idx,
                   samplesheet, verification_body):
    inventory = [
        ("genotypes/study_eur.bed", "study stratum 1 (EUR-drawn) genotypes"),
        ("genotypes/study_eur.bim", "study stratum 1 variants"),
        ("genotypes/study_eur.fam", "study stratum 1 samples"),
        ("genotypes/study_afr.bed", "study stratum 2 (AFR-drawn) genotypes"),
        ("genotypes/study_afr.bim", "study stratum 2 variants"),
        ("genotypes/study_afr.fam", "study stratum 2 samples"),
        ("reference/ref_panel.bed", "mini reference panel genotypes (autosomes only)"),
        ("reference/ref_panel.bim", "reference panel variants"),
        ("reference/ref_panel.fam", "reference panel samples"),
        ("reference/ref_pop.tsv", "`FID IID pop super_pop` for the reference panel"),
        ("samplesheet.csv", "pipeline input samplesheet"),
        ("highld/high_ld_regions_b38.bed", "GRCh38 long-range-LD intervals"),
    ]
    for pop in SUPER_POPS:
        inventory.append(("reference/keep/%s.keep" % pop,
                          "`FID IID` for the %s reference samples" % pop))
    inventory.sort()
    sizes = {}
    for rel, _desc in inventory:
        path = out / rel
        sizes[rel] = path.stat().st_size if path.exists() else 0

    per_chrom = []
    for c in range(1, 23):
        per_chrom.append((str(c), sum(1 for v in variants if v.chrom == str(c))))

    sample_summary = []
    for stratum in STRATA:
        samples, _geno = stratum_data[stratum["id"]]
        sample_summary.append((
            "genotypes/%s.fam" % stratum["id"], len(samples),
            sum(1 for s in samples if s.fam_sex == 1),
            sum(1 for s in samples if s.fam_sex == 2)))
    sample_summary.append((
        "reference/ref_panel.fam", len(ref_samples),
        sum(1 for s in ref_samples if s.fam_sex == 1),
        sum(1 for s in ref_samples if s.fam_sex == 2)))

    def ids(slots):
        """Slot IDs grouped by stratum: "`EUR_S020` + `EUR_S021`; `AFR_...`"."""
        return "; ".join(
            " + ".join("`%s_S%03d`" % (s["prefix"], slot) for slot in slots)
            for s in STRATA)

    def vids(vis):
        return ", ".join("`%s`" % variants[vi].vid for vi in vis)

    sample_defects = [
        ["a1", "Sex mismatch: `.fam` female, genotypes male",
         ids([SLOT_SEX_FAMFEMALE_GENOMALE]),
         "`.fam` sex = 2; chrX non-PAR calls are hemizygous (no hets) and all "
         "%d chrY variants are called" % sum(1 for v in variants if v.region == "Y"),
         "`--check-sex ycount` -> `PROBLEM`, X F ~ 1, YCOUNT high"],
        ["a2", "Sex mismatch: `.fam` male, genotypes female",
         ids([SLOT_SEX_FAMMALE_GENOFEMALE]),
         "`.fam` sex = 1; chrX non-PAR calls are diploid/heterozygous and every "
         "chrY call is missing",
         "`--check-sex ycount` -> `PROBLEM`, X F ~ 0, YCOUNT 0"],
        ["b1", "Duplicate / MZ pair",
         ids([SLOT_MZ_A, SLOT_MZ_B]),
         "the `S021` genotype vector is an exact copy of `S020` (all "
         "chromosomes)",
         "KING kinship ~ 0.5, IBS0 = 0 -> duplicate/MZ call"],
        ["b2", "First-degree pair (parent-offspring)",
         ids([SLOT_PARENT, SLOT_OFFSPRING]),
         "`S031` inherits one autosomal allele from `S030` and one from the "
         "population pool at every autosomal variant; no pedigree in the `.fam`",
         "KING kinship ~ 0.25 with IBS0 ~ 0 -> parent-offspring"],
        ["b3", "First-degree pair (full sibs)",
         ids([SLOT_SIB_A, SLOT_SIB_B]),
         "two virtual, non-genotyped parents; both sibs inherit one allele from "
         "each at every autosomal variant",
         "KING kinship ~ 0.25 with IBS0 > 0 -> full sibs"],
        ["c1", "High sample missingness",
         ids([40, 41]),
         "15% and 22% of all calls set missing at random",
         "`--missing` F_MISS > 0.10; removed by any `mind` <= 0.10"],
        ["c2", "Borderline sample missingness",
         ids([42]),
         "5% of all calls set missing at random",
         "F_MISS ~ 0.05: kept at `mind` 0.10, dropped at `mind` 0.02"],
        ["h", "Ancestry outliers (wrong super-population)",
         ", ".join("`%s_S%03d` (%s)" % (s["prefix"], slot, pop)
                   for s in STRATA
                   for slot, pop in sorted(SLOT_ANCESTRY_OUTLIER[s["id"]].items())),
         "genotypes drawn from another super-population's frequencies while the "
         "samplesheet labels the stratum EUR / AFR",
         "reference-projected PCA places them on the wrong cluster; ancestry "
         "assignment must not return the stratum label"],
    ]

    variant_defects = [
        ["d1", "High variant missingness (block of %d)" % len(defects["varmiss_high"]),
         vids(defects["varmiss_high"]),
         "chr3 block, %d-%d%% of samples set missing per variant (same "
         "variants in both strata)"
         % (int(VARMISS_HIGH_RANGE[0] * 100), int(VARMISS_HIGH_RANGE[1] * 100)),
         "`--missing` variant F_MISS > 0.10; removed by any `geno` <= 0.10"],
        ["d2", "Borderline variant missingness",
         vids(defects["varmiss_border"]),
         "chr4, a planted rate of %d-%d%% of samples per variant, rounded to a "
         "whole number of samples (hence the measured range below)"
         % (int(VARMISS_BORDER_RANGE[0] * 100), int(VARMISS_BORDER_RANGE[1] * 100)),
         "F_MISS ~ 0.05: kept at `geno` 0.10, dropped at `geno` 0.02"],
        ["e", "Gross HWE violation",
         vids(defects["hwe"]),
         "chr5, every sample forced heterozygous (observed het 1.00, expected 0.5)",
         "`--hardy` p ~ 0 -> excluded by any HWE threshold"],
        ["f", "MAF below 1%",
         vids(defects["rare"]),
         "chr6; every sample homozygous A2 except the first 1 (first three "
         "variants) or 2 (last three) of the fixed carriers %s -- "
         "MAF 1/250 = 0.40%% or 2/250 = 0.80%%" % ids(SLOT_RARE_CARRIERS),
         "`--freq` ALT_FREQS < 0.01 -> removed by `--maf 0.01`"],
        ["g", "Variants inside a listed high-LD interval",
         "see the high-LD table below",
         "a dense block of %d variants placed inside `hld_6_1_MHC` "
         "(`6:%d-%d`), plus whatever falls into the other intervals by chance"
         % (MHC_BLOCK_N, MHC_BLOCK[0], MHC_BLOCK[1]),
         "the high-LD exclusion step must drop them before LD pruning / PCA"],
    ]

    # ---- measured values -------------------------------------------------
    measured = []

    def per_stratum(fn):
        return [fn(*stratum_data[s["id"]]) for s in STRATA]

    def fmt(values, spec="%s"):
        return [spec % v for v in values]

    for slot in sorted(SLOT_MISS_SAMPLE):
        vals = per_stratum(lambda samples, geno, slot=slot: sample_f_miss(
            geno, next(i for i, s in enumerate(samples)
                       if s.iid.endswith("_S%03d" % slot))))
        measured.append(("F_MISS of `*_S%03d`" % slot,) + tuple(fmt(vals, "%.4f")))

    vals = per_stratum(lambda samples, geno: max(
        variant_stats(geno[vi])[1] for vi in defects["varmiss_high"]))
    measured.append(("max variant F_MISS in the chr3 block",) + tuple(fmt(vals, "%.4f")))
    vals = per_stratum(lambda samples, geno: min(
        variant_stats(geno[vi])[1] for vi in defects["varmiss_high"]))
    measured.append(("min variant F_MISS in the chr3 block",) + tuple(fmt(vals, "%.4f")))
    vals = per_stratum(lambda samples, geno: "%.4f-%.4f" % (
        min(variant_stats(geno[vi])[1] for vi in defects["varmiss_border"]),
        max(variant_stats(geno[vi])[1] for vi in defects["varmiss_border"])))
    measured.append(("borderline variant F_MISS range",) + tuple(fmt(vals)))
    vals = per_stratum(lambda samples, geno: "%.4f-%.4f" % (
        min(variant_stats(geno[vi])[0] for vi in defects["rare"]),
        max(variant_stats(geno[vi])[0] for vi in defects["rare"])))
    measured.append(("planted rare-variant MAF range",) + tuple(fmt(vals)))
    vals = per_stratum(lambda samples, geno: sum(
        1 for i, v in enumerate(variants)
        if v.region == "auto" and variant_stats(geno[i])[0] < 0.01))
    measured.append(("autosomal variants with MAF < 1% (planted + chance)",) + tuple(fmt(vals, "%d")))
    vals = per_stratum(lambda samples, geno: sum(
        1 for i, v in enumerate(variants)
        if v.region == "auto" and variant_stats(geno[i])[1] > 0.10))
    measured.append(("autosomal variants with F_MISS > 10%",) + tuple(fmt(vals, "%d")))
    vals = per_stratum(lambda samples, geno: sum(
        1 for i in range(len(samples)) if sample_f_miss(geno, i) > 0.10))
    measured.append(("samples with F_MISS > 10%",) + tuple(fmt(vals, "%d")))

    # ---- high-LD overlap -------------------------------------------------
    highld_overlap = []
    n_in_highld = 0
    for chrom, start, end, name in HIGH_LD_REGIONS_1BASED:
        count = sum(1 for v in variants
                    if v.chrom == str(chrom) and start <= v.pos <= end)
        n_in_highld += count
        if count:
            highld_overlap.append((name, "%d:%d-%d" % (chrom, start, end), count))
    n_in_mhc = sum(1 for v in variants
                   if v.chrom == "6" and MHC_BLOCK[0] <= v.pos <= MHC_BLOCK[1])

    verification = VERIFICATION.replace(
        "@@VERIFICATION_BODY@@",
        (verification_body or RECORDED_VERIFICATION).strip())

    return {
        "seed": args.seed,
        "inventory": inventory,
        "sizes": sizes,
        "n_var": len(variants),
        "n_var_auto": len(ref_variants),
        "n_par1": sum(1 for v in variants if v.region == "par1"),
        "n_nonpar": sum(1 for v in variants if v.region == "nonpar"),
        "n_par2": sum(1 for v in variants if v.region == "par2"),
        "n_y": sum(1 for v in variants if v.region == "Y"),
        "per_chrom": per_chrom,
        "sample_summary": sample_summary,
        "samplesheet": samplesheet,
        "n_flipped": len(flip_idx),
        "n_ambiguous": len(ambiguous_idx),
        "sample_defects": sample_defects,
        "variant_defects": variant_defects,
        "measured": measured,
        "highld_overlap": highld_overlap,
        "n_in_highld": n_in_highld,
        "n_in_mhc": n_in_mhc,
        "verification": verification,
        "id_lists": [
            ("Strand-flipped variants (study alleles complemented)",
             [variants[vi].vid for vi in flip_idx]),
            ("Strand-ambiguous variants (`A/T` or `C/G`, never flipped)",
             [variants[vi].vid for vi in ambiguous_idx]),
            ("Variants inside `hld_6_1_MHC`",
             [v.vid for v in variants
              if v.chrom == "6" and 25391794 <= v.pos <= 33424245]),
        ],
    }


if __name__ == "__main__":
    main()
