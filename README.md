# ![nfcore/test-datasets](docs/images/test-datasets_logo.png)

Test data to be used for automated testing with the nf-core pipelines

## Introduction

This is the `gwasqc` example-data branch. It holds the test fixtures for the
[gwasqc](https://github.com/lyh970817/gwasqc) genotype-QC pipeline, which ports
the GLAD/EDGI/NBR PLINK-format QC chain to nf-core-style Nextflow.

Like the other per-pipeline branches in this repository, this branch carries
only its own pipeline's data; it is not merged into `master` or `modules`.

## Git clone the gwasqc pipeline test data

```bash
# get a local copy of just this branch's data
git clone -b gwasqc --single-branch https://github.com/nf-core/test-datasets.git

# to update the data: fork first, then clone from your fork
git clone -b gwasqc --single-branch https://github.com/USERNAME/test-datasets.git
```

## Contents

```
data/gwasqc/
  genotypes/<stratum>.{bed,bim,fam}       two study strata (study_eur, study_afr)
  reference/ref_panel.{bed,bim,fam}       mini five-super-population reference panel
  reference/ref_pop.tsv                   FID IID pop super_pop
  reference/keep/<SUPERPOP>.keep          FID IID per super-population
  highld/high_ld_regions_b38.bed          GRCh38 long-range-LD intervals
  samplesheet.csv                         pipeline input (id,ancestry,bed,bim,fam)
  README.md                               every planted defect, its IDs and its
                                          expected QC outcome, plus the
                                          verification commands and results
scripts/gwasqc/generate_fixtures.py       deterministic generator (fixed seed,
                                          Python standard library only)
```

**The data is synthetic.** No real human genotypes are present. Everything under
`data/gwasqc/` is reproduced from scratch by

```bash
python3 scripts/gwasqc/generate_fixtures.py
```

Start with [`data/gwasqc/README.md`](data/gwasqc/README.md): it is the contract
for what these fixtures contain and what a correct QC pipeline is expected to
find in them.

## Documentation

nf-core/test-datasets comes with documentation in the `docs/` directory:

01. [Add a new  test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/ADD_NEW_DATA.md)
02. [Use an existing test dataset](https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md)

## Support

For further information or help, don't hesitate to get in touch on our [Slack organisation](https://nf-co.re/join/slack) (a tool for instant messaging).
