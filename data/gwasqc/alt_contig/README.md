# Alt- and random-contig copies of a stratum and the reference panel

`genotypes/study_eur` and `reference/ref_panel` from the parent directory, each with eight variants appended
on GRCh38 alt and random contigs, to exercise the exclusion of variants whose chromosome code is not one of
1-22, X, Y, XY and MT. Per-ancestry Axiom array filesets on GRCh38 carry such variants in their `.bim`.

| Path | Contents |
| --- | --- |
| `study_eur_alt.{bed,bim,fam}` | `study_eur` (3962 variants, 125 samples) followed by the 8 contig variants: 3970 variants |
| `ref_panel_alt.{bed,bim,fam}` | `ref_panel` (3632 variants, 250 samples) followed by the same 8 contig variants: 3640 variants |

The first 3962 (respectively 3632) variants of each fileset, their `.bed` blocks and the `.fam` are those of
the original unchanged, byte for byte. Each contig variant is a copy of one autosomal variant of the same
fileset, genotypes and alleles unchanged, under a new ID, chromosome code and position:

| Copied from | Chromosome code | ID | Position |
| --- | --- | --- | --- |
| `rs1000514` | `1_KI270706v1_random` | `AX-900000001` | 50001 |
| `rs1802793` | `7_KI270803v1_alt` | `AX-900000002` | 120001 |
| `rs1802804` | `7_KI270803v1_alt` | `AX-900000003` | 240001 |
| `rs1901097` | `8_KI270821v1_alt` | `AX-900000004` | 30001 |
| `rs2592560` | `17_KI270909v1_alt` | `AX-900000005` | 90001 |
| `rs2696237` | `19_KI270938v1_alt` | `AX-900000006` | 60001 |
| `rs2800894` | `22_KI270879v1_alt` | `AX-900000007` | 150001 |
| `rs2800909` | `22_KI270928v1_alt` | `AX-900000008` | 80001 |

Each source variant passes `--maf 0.01`, `--geno 0.05` and `--hwe 1e-10` in `study_eur` and is neither
strand-flipped nor strand-ambiguous, so a copy that reached the variant filters would survive them all: only a
contig exclusion removes it. The study and the reference carry the copies under the same IDs, positions and
alleles.

Reproduce with plink2 2.0.0a.6.9 and coreutils from the parent directory (a `.bed` is variant-major, so the
contig variants' blocks are the restricted fileset's blocks after its 3-byte header):

```bash
cat > contigs.tsv <<'EOF'
rs1000514	1_KI270706v1_random	AX-900000001	50001
rs1802793	7_KI270803v1_alt	AX-900000002	120001
rs1802804	7_KI270803v1_alt	AX-900000003	240001
rs1901097	8_KI270821v1_alt	AX-900000004	30001
rs2592560	17_KI270909v1_alt	AX-900000005	90001
rs2696237	19_KI270938v1_alt	AX-900000006	60001
rs2800894	22_KI270879v1_alt	AX-900000007	150001
rs2800909	22_KI270928v1_alt	AX-900000008	80001
EOF
cut -f1 contigs.tsv > sources.txt
for pair in genotypes/study_eur:study_eur_alt reference/ref_panel:ref_panel_alt; do
  in=${pair%%:*}; out=alt_contig/${pair##*:}
  plink2 --bfile "$in" --extract sources.txt --make-bed --out sources
  awk 'BEGIN { OFS = "\t" } NR == FNR { chrom[$1] = $2; id[$1] = $3; pos[$1] = $4; next } { print chrom[$2], id[$2], $3, pos[$2], $5, $6 }' contigs.tsv sources.bim > contigs.bim
  cat "$in.bim" contigs.bim > "$out.bim"
  { cat "$in.bed"; tail -c +4 sources.bed; } > "$out.bed"
  cp "$in.fam" "$out.fam"
done
```
