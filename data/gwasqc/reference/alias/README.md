# Alias-named, split reference panel

`ref_panel` from the parent directory, split into two filesets and renamed, to exercise position-based
harmonisation of a reference panel whose variant IDs differ from the study's.

| Path | Contents |
| --- | --- |
| `ref_alias_a.{bed,bim,fam}` | `ref_panel` chromosomes 1-11 (2464 variants), bare chromosome codes |
| `ref_alias_b.{bed,bim,fam}` | `ref_panel` chromosomes 12-22 (1168 variants), chromosome codes written `chr12`-`chr22` |

Genotypes, samples, positions and alleles are those of `ref_panel` unchanged. In each `.bim`, every seventh
variant (rows 7, 14, ...) carries the ID `AX-<position, zero-padded to 9 digits>` instead of its rsID: 352
variants in `ref_alias_a`, 166 in `ref_alias_b`, 518 in all. Every ID is unique across the two files. Matched to
a study stratum by ID alone these 518 variants are lost; matched by chromosome, position and alleles every
variant of `ref_panel` is recovered, and the merged panel, renamed to the study's IDs, is identical to
`ref_panel` restricted to the same variants.

Reproduce with plink2 2.0.0a.6.9 and awk from the parent directory:

```bash
plink2 --bfile ref_panel --chr 1-11 --make-bed --out a
plink2 --bfile ref_panel --chr 12-22 --make-bed --out b
cp a.bed alias/ref_alias_a.bed; cp b.bed alias/ref_alias_b.bed
cp ref_panel.fam alias/ref_alias_a.fam; cp ref_panel.fam alias/ref_alias_b.fam
awk 'BEGIN{OFS="\t"} {if (NR % 7 == 0) $2 = "AX-" sprintf("%09d", $4); print}' a.bim > alias/ref_alias_a.bim
awk 'BEGIN{OFS="\t"} {if (NR % 7 == 0) $2 = "AX-" sprintf("%09d", $4); $1 = "chr" $1; print}' b.bim > alias/ref_alias_b.bim
```
