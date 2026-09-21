# gwasqc test fixtures

Synthetic genotype fixtures for the `gwasqc` genotype-QC pipeline (<https://github.com/lyh970817/gwasqc>).
**No real human data.** Every genotype here is simulated.

Everything in this directory is reproduced from scratch by

```bash
python3 scripts/gwasqc/generate_fixtures.py   # seed 20260920, stdlib only
```

The generator is deterministic: a re-run with the same seed reproduces byte-identical files. Changing the script changes the bytes, so treat the committed data as the artefact and the script as its provenance.

## Files

| Path | Size | Contents |
| --- | --- | --- |
| `genotypes/study_afr.bed` | 123.8 KB | study stratum 2 (AFR-drawn, 125 samples) genotypes |
| `genotypes/study_afr.bim` | 106.8 KB | study stratum 2 variants |
| `genotypes/study_afr.fam` | 3.3 KB | study stratum 2 samples |
| `genotypes/study_csa.bed` | 11.6 KB | study stratum 3 (CSA-drawn, 12 samples) genotypes |
| `genotypes/study_csa.bim` | 106.8 KB | study stratum 3 variants |
| `genotypes/study_csa.fam` | 324 B | study stratum 3 samples |
| `genotypes/study_eur.bed` | 123.8 KB | study stratum 1 (EUR-drawn, 125 samples) genotypes |
| `genotypes/study_eur.bim` | 106.8 KB | study stratum 1 variants |
| `genotypes/study_eur.fam` | 3.3 KB | study stratum 1 samples |
| `highld/high_ld_regions_b38.bed` | 715 B | GRCh38 long-range-LD intervals |
| `reference/keep/AFR.keep` | 800 B | `FID IID` for the AFR reference samples |
| `reference/keep/AMR.keep` | 800 B | `FID IID` for the AMR reference samples |
| `reference/keep/CSA.keep` | 800 B | `FID IID` for the CSA reference samples |
| `reference/keep/EAS.keep` | 800 B | `FID IID` for the EAS reference samples |
| `reference/keep/EUR.keep` | 800 B | `FID IID` for the EUR reference samples |
| `reference/ref_panel.bed` | 223.5 KB | mini reference panel genotypes (autosomes only) |
| `reference/ref_panel.bim` | 98.1 KB | reference panel variants |
| `reference/ref_panel.fam` | 6.1 KB | reference panel samples |
| `reference/ref_pop.tsv` | 5.9 KB | `FID IID pop super_pop` for the reference panel |
| `samplesheet.csv` | 282 B | pipeline input samplesheet |

Total, excluding this README: **925.1 KB**.

## Conventions

* **Build: GRCh38 / b38.** Positions are GRCh38 primary-assembly coordinates; chrX PAR boundaries match plink1.9 `--split-x b38` and plink2 `--split-par b38` (PAR1 `X:10001-2781479`, PAR2 `X:155701383-156030895`).
* **Chromosome codes** in `.bim` are bare PLINK codes without a `chr` prefix: `1`-`22`, `X`, `Y` (plus `XY` after a PAR split). `high_ld_regions_b38.bed` uses the same codes so it can be joined straight onto a `.bim`/`.pvar`.
* `.bim` is tab-delimited, `.fam` space-delimited (what PLINK itself writes). `.bim` column 5 is A1 (the counted allele), column 6 is A2.
* `.fam`/`ref_panel.fam` phenotype column is `-9` throughout: this is a QC fixture, not an association fixture.
* `FID == IID` for every sample, and the `.fam` PAT/MAT columns are `0`. The related pairs below carry **no pedigree**: relatedness has to be inferred from genotypes, which is the point.
* `high_ld_regions_b38.bed` is a real BED: **0-based, half-open**. Its interval set is identical to the legacy asset `GLAD_EDGI_NBR_qc_workflows/full_pipe/dependencies/highLDregions4bim_b38.awk` (which states the same intervals 1-based inclusive).

## Composition

### Variants (3962 in the study filesets, 3632 in the reference)

| Region | Variants |
| --- | --- |
| autosomes chr1-22 | 3632 |
| chrX PAR1 (`X:10001-2781479`) | 30 |
| chrX non-PAR | 250 |
| chrX PAR2 (`X:155701383-156030895`) | 20 |
| chrY (male-specific region) | 30 |

Per-autosome counts: chr1 312, chr2 303, chr3 248, chr4 238, chr5 227, chr6 244, chr7 200, chr8 182, chr9 173, chr10 168, chr11 169, chr12 167, chr13 143, chr14 134, chr15 128, chr16 113, chr17 104, chr18 101, chr19 73, chr20 81, chr21 60, chr22 64.

Variant IDs are rsID-style and strictly increasing with genomic order. The 50 chrX PAR variants exist so that the PAR split has work to do: both plink1.9 `--split-x b38` and plink2 `--split-par b38` move exactly these 50 variants to `XY`/`25`, leaving the 250 non-PAR chrX variants that the sex check uses. plink2 has **no** `--split-x`; there the flag is spelled `--split-par`.

### Samples

| Fileset | Samples | `.fam` sex 1 / 2 |
| --- | --- | --- |
| `genotypes/study_eur.fam` | 125 | 64 / 61 |
| `genotypes/study_afr.fam` | 125 | 64 / 61 |
| `genotypes/study_csa.fam` | 12 | 5 / 7 |
| `reference/ref_panel.fam` | 250 | 130 / 120 |

Study sample IDs are `<PREFIX>_S<nnn>` with `PREFIX` = `EUR` and `nnn` = 001-125 in `study_eur`, `PREFIX` = `AFR` and `nnn` = 001-125 in `study_afr`, `PREFIX` = `CSA` and `nnn` = 001-012 in `study_csa`; reference sample IDs are `<SUBPOP>_<nnn>`. The numeric part of a study ID is stable across regenerations, so the defect IDs below never move.

### `samplesheet.csv`

```csv
id,ancestry,bed,bim,fam
study_eur,EUR,genotypes/study_eur.bed,genotypes/study_eur.bim,genotypes/study_eur.fam
study_afr,AFR,genotypes/study_afr.bed,genotypes/study_afr.bim,genotypes/study_afr.fam
study_csa,CSA,genotypes/study_csa.bed,genotypes/study_csa.bim,genotypes/study_csa.fam
```

Header is exactly `id,ancestry,bed,bim,fam`; the three path columns are **relative to this directory** (`data/gwasqc/`).

## Ancestry design

Allele frequencies come from a two-level Balding-Nichols model: one ancestral frequency per variant drawn from U(0.10, 0.50), differentiated across the five super-populations with **Fst = 0.12**, then across two sub-populations per super-population with Fst = 0.01. Frequencies are clamped to [0.03, 0.97]. Genotypes are binomial draws, so every variant is HWE-conforming unless a defect was planted into it.

| Super-population | Sub-populations | Reference samples |
| --- | --- | --- |
| AFR | YRI, ESN | 50 |
| AMR | MXL, PEL | 50 |
| CSA | GIH, PJL | 50 |
| EAS | CHB, JPT | 50 |
| EUR | GBR, TSI | 50 |

`reference/ref_pop.tsv` is tab-separated with the header `FID IID pop super_pop`; `reference/keep/<SUPERPOP>.keep` lists `FID IID` (space-separated, no header) for that super-population.

Study strata are drawn from their matching super-population frequencies (`study_eur` from EUR, `study_afr` from AFR, `study_csa` from CSA), so a reference-projected PCA should place them on top of the matching reference cluster -- except for the planted ancestry outliers.

**Reference/study harmonisation.** Reference and study filesets share variant ID, position and alleles on all 3632 autosomal variants, with two documented exceptions that exist to exercise harmonisation:

* **54 strand-flipped variants** (1.5% of autosomal variants): the study `.bim` carries the complement of the reference alleles (e.g. reference `A G` -> study `T C`) while the genotype dosages are unchanged, which is exactly what a strand flip looks like. All of them have unambiguous (non-`A/T`, non-`C/G`) allele pairs, so the flip is resolvable. Every stratum flips the same variant set.
* **40 strand-ambiguous variants** carry `A/T` or `C/G` alleles in both the reference and the study sets (never flipped). These are the ones an IUPAC/ambiguity check must refuse to resolve.

Full ID lists: see `flipped variants` and `ambiguous variants` at the bottom of this file.

## Planted defects

Every ID below is deterministic: the slot number in a study sample ID encodes its role. The two 125-sample strata carry the same slots; the 12-sample stratum carries only the defects listed against its IDs. Sample-level defect IDs are listed per stratum, for the strata that carry the defect; variant-level defects hit the same variant IDs in every stratum.

### Sample-level

| # | Defect | IDs | How it was planted | Expected QC outcome |
| --- | --- | --- | --- | --- |
| a1 | Sex mismatch: `.fam` female, genotypes male | `EUR_S010`; `AFR_S010`; `CSA_S010` | `.fam` sex = 2; chrX non-PAR calls are hemizygous (no hets) and all 30 chrY variants are called | `--check-sex ycount` -> `PROBLEM`, X F ~ 1, YCOUNT high |
| a2 | Sex mismatch: `.fam` male, genotypes female | `EUR_S011`; `AFR_S011` | `.fam` sex = 1; chrX non-PAR calls are diploid/heterozygous and every chrY call is missing | `--check-sex ycount` -> `PROBLEM`, X F ~ 0, YCOUNT 0 |
| b1 | Duplicate / MZ pair | `EUR_S020` + `EUR_S021`; `AFR_S020` + `AFR_S021`; `CSA_S005` + `CSA_S006` | the second sample's genotype vector is an exact copy of the first's (all chromosomes) | KING kinship ~ 0.5, IBS0 = 0 -> duplicate/MZ call |
| b2 | First-degree pair (parent-offspring) | `EUR_S030` + `EUR_S031`; `AFR_S030` + `AFR_S031` | the second sample inherits one autosomal allele from the first and one from the population pool at every autosomal variant; no pedigree in the `.fam` | KING kinship ~ 0.25 with IBS0 ~ 0 -> parent-offspring |
| b3 | First-degree pair (full sibs) | `EUR_S032` + `EUR_S033`; `AFR_S032` + `AFR_S033` | two virtual, non-genotyped parents; both sibs inherit one allele from each at every autosomal variant | KING kinship ~ 0.25 with IBS0 > 0 -> full sibs |
| c1 | High sample missingness | `EUR_S040` + `EUR_S041`; `AFR_S040` + `AFR_S041` | 15% and 22% of all calls set missing at random | `--missing` F_MISS > 0.10; removed by any `mind` <= 0.10 |
| c2 | Borderline sample missingness | `EUR_S042`; `AFR_S042` | 5% of all calls set missing at random | F_MISS ~ 0.05: kept at `mind` 0.10, dropped at `mind` 0.02 |
| c3 | Sparse sample: one missing call at a handful of variants | `CSA_S003` | exactly one genotype set missing at 5 autosomal variants outside every planted block (`rs1151545`, `rs1311999`, `rs1911972`, `rs2167038`, `rs2770893`), which adds 5/3962 to the sample's own F_MISS (the measured value below includes the calls the chr3/chr4 blocks assign to it at random) | the sample survives any `mind`; each of the 5 variants has F_MISS 1/n = 0.0833 in a 12-sample stratum, above `geno` 0.05, so they fall at geno although only one call is missing: the 1/n granularity of missingness at small n |
| h | Ancestry outliers (wrong super-population) | `EUR_S050` (EAS), `EUR_S051` (AFR), `EUR_S052` (AMR), `AFR_S050` (EUR), `AFR_S051` (CSA), `AFR_S052` (EAS) | genotypes drawn from another super-population's frequencies while the samplesheet labels the stratum EUR / AFR | reference-projected PCA places them on the wrong cluster; ancestry assignment must not return the stratum label |

### Variant-level

| # | Defect | IDs | How it was planted | Expected QC outcome |
| --- | --- | --- | --- | --- |
| d1 | High variant missingness (block of 20) | `rs1395803`, `rs1396390`, `rs1396651`, `rs1397483`, `rs1397593`, `rs1398175`, `rs1398777`, `rs1399746`, `rs1400520`, `rs1401325`, `rs1401932`, `rs1402148`, `rs1403097`, `rs1403250`, `rs1403949`, `rs1404238`, `rs1405053`, `rs1405625`, `rs1405639`, `rs1406332` | chr3 block, 12-25% of samples set missing per variant (same variants in both strata) | `--missing` variant F_MISS > 0.10; removed by any `geno` <= 0.10 |
| d2 | Borderline variant missingness | `rs1498627`, `rs1499588`, `rs1500478`, `rs1501405`, `rs1502040` | chr4, a planted rate of 4-6% of samples per variant, rounded to a whole number of samples (hence the measured range below) | F_MISS ~ 0.05: kept at `geno` 0.10, dropped at `geno` 0.02 |
| e | Gross HWE violation | `rs1667531`, `rs1668208`, `rs1668410`, `rs1668428`, `rs1668730` | chr5, every sample forced heterozygous (observed het 1.00, expected 0.5) | `--hardy` p ~ 0 -> excluded by any HWE threshold |
| f | MAF below 1% | `rs1746871`, `rs1746964`, `rs1747289`, `rs1748137`, `rs1748739`, `rs1749463` | chr6; every sample homozygous A2 except the first 1 (first three variants) or 2 (last three) of the fixed carriers `EUR_S060` + `EUR_S061`; `AFR_S060` + `AFR_S061` -- MAF 1/250 = 0.40% or 2/250 = 0.80%; no carrier at all in `study_csa`, where the six variants are monomorphic (MAF 0) | `--freq` ALT_FREQS < 0.01 -> removed by `--maf 0.01` |
| g | Variants inside a listed high-LD interval | see the high-LD table below | a dense block of 30 variants placed inside `hld_6_1_MHC` (`6:25400000-33400000`), plus whatever falls into the other intervals by chance | the high-LD exclusion step must drop them before LD pruning / PCA |

### Small stratum

`study_csa` has 12 samples, so every per-sample and per-variant rate is a multiple of 1/12 = 0.0833. Three consequences, all measured in the committed data and none of them defects of the fixture:

* **One missing call is 8.3% of the samples**, above `geno` 0.05: the 5 sparse-sample variants (defect c3) and every variant of the chr4 borderline block (defect d2, whose planted 4-6% rounds to one missing call at n = 12) fall at `geno` here, while the borderline block survives in the 125-sample strata. A per-variant rate is a fraction of the samples, so this granularity is a property of the sample count; a per-sample rate (`mind`) is a fraction of the variants and is unaffected by it.
* **The HWE plants (defect e) are not detectable at n = 12.** With every one of 12 samples heterozygous the exact-test p-value is of the order of 1e-3 (the all-heterozygous table itself has probability 1.5e-03 under HWE; plink2's value is in the verification below), nowhere near a 1e-10 threshold; `--hwe 1e-10` therefore excludes nothing in this stratum.
* **MAF below 1% means monomorphic at n = 12** (the smallest non-zero MAF is 1/24 = 0.0417), so the planted rare variants (defect f) are written with no carrier and a number of unplanted variants are monomorphic by sampling chance; the measured count below is the real total that `--maf 0.01` removes.

Its 12 samples also make it the stratum on which the `gwasqc` small-n route (reference-anchored sample QC) is exercised: LD pruning on so few samples is the failure the route exists to avoid.

### Measured values in the committed data

| Statistic | `study_eur` | `study_afr` | `study_csa` |
| --- | --- | --- | --- |
| F_MISS of `*_S040` | 0.1499 | 0.1499 | - |
| F_MISS of `*_S041` | 0.2201 | 0.2201 | - |
| F_MISS of `*_S042` | 0.0500 | 0.0500 | - |
| F_MISS of the sparse sample `*_S003` | - | - | 0.00278 |
| max variant F_MISS in the chr3 block | 0.2560 | 0.2480 | 0.2500 |
| min variant F_MISS in the chr3 block | 0.1600 | 0.1280 | 0.0833 |
| borderline variant F_MISS range | 0.0480-0.0640 | 0.0480-0.0640 | 0.0833-0.0833 |
| planted rare-variant MAF range | 0.0040-0.0081 | 0.0040-0.0081 | 0.0000-0.0000 |
| autosomal variants with MAF < 1% (planted + chance) | 8 | 7 | 236 |
| autosomal variants with F_MISS > 10% | 20 | 20 | 18 |
| autosomal variants with F_MISS > 5% | 24 | 23 | 30 |
| samples with F_MISS > 10% | 2 | 2 | 0 |

## High-LD regions

`highld/high_ld_regions_b38.bed` lists 24 long-range-LD intervals (GRCh38), tab-separated, 4 columns `chrom start end name`, 0-based half-open.

**202 of the 3632 autosomal variants fall inside a listed interval**, of which 38 sit in the deliberately dense block inside `hld_6_1_MHC` (`6:25400000-33400000`). Overlap counts per interval that catches anything:

| Interval | Region | Variants |
| --- | --- | --- |
| `hld_1_1` | `1:47822309-51822307` | 4 |
| `hld_2_1` | `2:85861220-100425020` | 15 |
| `hld_2_2_LCT` | `2:133908698-137408698` | 5 |
| `hld_2_3` | `2:182309768-189309768` | 7 |
| `hld_3_1` | `3:47483507-49987563` | 6 |
| `hld_3_2` | `3:83368160-86868160` | 2 |
| `hld_3_3` | `3:88868161-96298466` | 6 |
| `hld_5_1` | `5:44464142-51168409` | 13 |
| `hld_5_2` | `5:98636397-101136397` | 4 |
| `hld_5_3` | `5:129636409-132636409` | 5 |
| `hld_5_4` | `5:136136413-139136412` | 3 |
| `hld_6_1_MHC` | `6:25391794-33424245` | 38 |
| `hld_6_2` | `6:57027244-63232136` | 11 |
| `hld_6_3` | `6:139637171-142137170` | 4 |
| `hld_7_1` | `7:55158099-67090863` | 9 |
| `hld_8_1_INV8p23` | `8:8105069-12105082` | 3 |
| `hld_8_2` | `8:43025701-48924888` | 11 |
| `hld_8_3` | `8:110918596-113918595` | 7 |
| `hld_10_1` | `10:36671067-43184546` | 14 |
| `hld_11_1` | `11:46021874-57475951` | 15 |
| `hld_11_2` | `11:88127185-91127184` | 5 |
| `hld_12_1` | `12:32955800-41319931` | 8 |
| `hld_12_2` | `12:110599476-113099475` | 2 |
| `hld_20_1` | `20:33948534-36438183` | 5 |

## Verification

Verified with the pinned containers

* `quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0` (PLINK v2.0.0-a.6.9LM, 29 Jan 2025)
* `quay.io/biocontainers/plink:1.90b6.21--h779adbc_1` (PLINK v1.90b6.21, 19 Oct 2020)

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

### 11. `study_csa` -- the 12-sample stratum

Verified separately, with
`community.wave.seqera.io/library/plink2:2.0.0a.6.9--e6710830a4b7f0c6` (the
same PLINK v2.0.0-a.6.9LM build, the image the pipeline pins); plink1.9 was
not run on this stratum.

```
plink2 --bfile genotypes/study_csa --freq --missing --hardy --out out_csa
```

`12 samples (7 females, 5 males)`, `3962 variants`.  The six planted chr6
rare variants have `ALT_FREQS 0`, `OBS_CT 24`; 236 autosomal variants are
monomorphic in this stratum (the 6 plants plus 230 by sampling chance) and no
other autosomal variant has MAF < 0.01, so `--maf 0.01` removes 236.  The
five sparse-sample variants (`rs1151545`, `rs1311999`, `rs1911972`,
`rs2167038`, `rs2770893`) each have `MISSING_CT 1`, `F_MISS 0.0833`, and
`--export A` shows the missing call is `CSA_S003`'s at all five; `CSA_S003`
itself has `F_MISS 0.00278` (11 of 3962: the 5 planted calls plus 6 that the
chr3/chr4 blocks assign to it at random).  30 autosomal variants have
`F_MISS > 0.05`: the 20 chr3 block variants (0.0833-0.25), the 5 chr4
borderline variants (each exactly one missing call, 0.0833) and the 5
sparse-sample variants.  The five chr5 HWE plants have `O(HET_A1) 1`,
`E(HET_A1) 0.5`, `P 0.0018564`, and the smallest autosomal `P` in the stratum
is 6.7e-4 (`rs2009918`), so `--hwe 1e-10` excludes nothing here.

**Small-sample guards.** plink2 2.0.0a.6.9 refuses to impute allele
frequencies or LD from fewer than 50 samples: `--check-sex`, `--het` and
`--pca` stop with `Error: This run requires decent allele frequencies, but
they aren't being loaded with --read-freq, and less than 50 samples are
available to impute them from.`, and `--indep-pairwise` stops with `Error:
This run estimates linkage disequilibrium between variants, but there are less
than 50 samples to estimate from.`  The overrides are `--bad-freqs` for the
first three and `--bad-ld` for pruning; a `--read-freq <file>` also lifts the
frequency guard for the whole run, even when the file covers only the
autosomes (chrX frequencies are then imputed from the 12 samples without
further notice).

Sex check:

```
plink2 --bfile genotypes/study_csa --split-par b38 --make-bed --out p2sx_csa
plink2 --bfile p2sx_csa --bad-freqs --check-sex max-female-xf=0.2 min-male-xf=0.8 \
       max-female-ycount=0 min-male-ycount=1 cols=+ycount --out p2cs_csa
```

`--split-par: 50 chromosome codes changed`, then `--check-sex: 250 chrX
variants and 30 variants scanned, 1 problem detected`:

```
FID       IID       PEDSEX SNPSEX STATUS  F  YCOUNT YRATE
CSA_S010  CSA_S010  2      1      PROBLEM 1  30     1
```

The females' F is 0.014 to 0.126 (below `max-female-xf` 0.2) and every
genotype male has F 1 and YCOUNT 30.  `--read-freq` of the stratum's own
`--freq` output, or of the autosomal CSA-reference frequencies, gives the same
12 verdicts (F differs in the sixth decimal).  22 monomorphic chrX variants
are skipped with a warning.

Relatedness:

```
plink2 --bfile genotypes/study_csa --chr 1-22 --make-king-table --king-table-filter 0.08 --out king_csa
```

`1 relationship reported (65 filtered out)`:

```
IID1      IID2      NSNP  HETHET    IBS0  KINSHIP
CSA_S006  CSA_S005  3623  0.360199  0     0.5
```

No other pair reaches 0.08.  KING needs no allele frequencies and runs bare
at n = 12.

Heterozygosity: `--het` bare is refused (above).  `--bad-freqs --het` gives
mean F -0.046 over the 12 samples (the in-sample frequency estimate is biased
at n = 12); `--het --read-freq` of the 50 CSA reference samples' `--freq`
over the 3578 shared same-strand variants gives mean F -0.0006, and the
reference samples themselves on that basis -0.0078, so a matched reference
frequency file puts study and reference on one scale.

LD pruning: `--bad-ld --indep-pairwise 1500 150 0.2 --chr 1-22` removes 3459
of 3632 autosomal variants and keeps 173 (4.8%), on data with no LD structure
at all: pairwise r^2 over 12 samples is inflated enough to prune almost
everything.  This is the small-n pruning collapse the stratum exists to
demonstrate; the 125-sample strata keep 3617 of 3632 with the same settings.

Ancestry, by the pipeline's own projection recipe: reference
`--maf 0.01 --freq counts --pca 6 allele-wts vcols=chrom,ref,alt` over the
3578 autosomal variants whose alleles agree between study and reference (the
54 documented flips excluded), then
`--score <weights> 2 5 header-read no-mean-imputation variance-standardize --score-col-nums 6-11 --read-freq <counts>`
for the reference and for `study_csa` alike.  Nearest-centroid assignment on
PC1-PC4: reference self-assignment 250/250, and all 12 `study_csa` samples
nearest the CSA centroid (distance 0.011-0.033).

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
* **`study_csa` was checked with plink2 only** (section 11); the plink1.9
  route of sections 2 and 8 was not repeated on it.


## Simplifications and caveats

* The reference panel is **autosomes only** (3632 variants). The ancestry PCA only uses autosomes; nothing here can test a reference chrX path.
* `AMR` and `CSA` are modelled as plain Balding-Nichols populations, **not** as admixed ones. Real AMR/CSA samples are admixed; a classifier tuned on this fixture will look better than it is on real data.
* There is **no LD structure**: variants are drawn independently. A high-LD region here is high-LD only by coordinate, not by correlation, so `--indep-pairwise` will prune almost nothing and the high-LD exclusion asset is exercised by coordinate overlap alone.
* Defects are planted in the **study strata only**; the reference panel is clean (no missingness, no HWE violation, no relatives, no sex mismatch).
* Relatedness is planted on the **autosomes**; the derived samples' chrX non-PAR and chrY genotypes are drawn independently.
* Hemizygous calls (chrX non-PAR and chrY, for samples whose *genotypes* are male) are written as homozygotes, which is PLINK's haploid convention. Two PLINK warnings are therefore **expected and not fixture defects**: before the PAR split every male's diploid call at the 50 PAR variants counts as `het. haploid` (~1.2k calls), and after the split the ~80 that remain are defect a2's heterozygous chrX calls; `Nonmissing nonmale Y chromosome genotype(s) present` is defect a1.
* Planted rare variants are set to a fixed carrier count rather than sampled, so their MAF is exact. A handful of *unplanted* variants also fall below MAF 1% by sampling chance; the measured counts above give the real totals.
* The 12-sample stratum is small by design and behaves as small strata do: see *Small stratum* above for what its size does to the missingness, HWE and MAF filters.

## Appendix: variant ID lists

<details><summary>Strand-flipped variants (study alleles complemented) (54)</summary>

```
rs1019157 rs1084781 rs1090252 rs1097655 rs1145782 rs1159607 rs1163718 rs1273151
rs1360139 rs1389123 rs1436911 rs1452101 rs1467337 rs1549826 rs1578801 rs1677459
rs1680161 rs1722029 rs1723746 rs1762883 rs1776206 rs1803389 rs1862912 rs1916892
rs1941077 rs1988230 rs1996775 rs2027882 rs2109777 rs2130274 rs2144516 rs2150265
rs2164106 rs2169905 rs2281679 rs2294093 rs2309898 rs2314061 rs2347169 rs2362593
rs2389678 rs2439958 rs2469140 rs2504154 rs2515327 rs2516136 rs2517896 rs2524188
rs2526561 rs2637925 rs2675388 rs2729284 rs2781271 rs2789194
```

</details>

<details><summary>Strand-ambiguous variants (`A/T` or `C/G`, never flipped) (40)</summary>

```
rs1008974 rs1111502 rs1117409 rs1212885 rs1239733 rs1301534 rs1323806 rs1325713
rs1462042 rs1515061 rs1532976 rs1540256 rs1617881 rs1630597 rs1648278 rs1658669
rs1676231 rs1678804 rs1685698 rs1719304 rs1801805 rs1810040 rs1880092 rs2032830
rs2036157 rs2115243 rs2140434 rs2173298 rs2239178 rs2293247 rs2339427 rs2417072
rs2489817 rs2586694 rs2670145 rs2679890 rs2693456 rs2696100 rs2779816 rs2814949
```

</details>

<details><summary>Variants inside `hld_6_1_MHC` (38)</summary>

```
rs1691803 rs1692225 rs1692986 rs1693180 rs1693565 rs1694305 rs1695168 rs1695930
rs1696394 rs1696687 rs1697113 rs1697536 rs1697583 rs1698450 rs1699148 rs1699176
rs1699973 rs1700763 rs1701059 rs1701798 rs1701848 rs1702222 rs1702918 rs1703899
rs1704251 rs1704269 rs1704945 rs1705262 rs1705437 rs1706350 rs1706917 rs1707913
rs1708405 rs1708542 rs1709245 rs1709713 rs1710071 rs1710602
```

</details>

