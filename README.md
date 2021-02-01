# PolygenicRSS


![DAG or simulation](wf.png)

# Configuration

first, download the repo do a directory of your choosing, then `cd PolygenicRSS/workflow`.  You will find an elaborate `config_base.yaml`.
You may find it useful to modify this file to run the pipeline outside of `cri`.


# Running the pipeline

## A primer on `snakemake`

You can read a lot more about snakemake [here](https://snakemake.readthedocs.io/en/stable/), but briefly, 
snakemake is a tool to create reproducible and scalable data analyses.  When in a directory with a file called 
`snakefile`, (like `PolygenicRSS/workflow`), and `snakemake` is on the path (via something like `spack load snakemake`), 
`snakemake -j 350 --profile ../config/gardner results/sim_ukb_ind/polym_3_10_0.5_10000_1.int.log"` will 
instruct snakemake to create the file `./results/sim_ukb_ind/polym_3_10_0.5_10000_1.int.log` using up to 350 simultaneous cluster jobs, with configuration about the cluster (and other things, in the directory `../config/gardner/`.  This will cause ld score regression to be run on the gwas for a trait with an `h2` of `0.5`, specifically on the `3`rd replicate of that gwas, with causal variants coming from the `polym` set of SNPs.  The GWAS will have 10000 individuals, and LD score regression will be run with LD scores generated with a 1cm window, and LD score regression itself will be run allowing the intercept term to vary.  By changing the file name that we ask snakemake to create, we can generate a different set of results.

## Snakemake rule syntax.

snakemake will go through all the "rules" it knows about to try to create the file.  It starts in `snakefile`, and in `snakefile` we see the line  `include: "ukb_snakefile"`.  Looking in the `ukb_snakefile` file we seethe following rule:

```snakemake
  rule baseline_ldsc:
      input:
          gwasf=config_d['GWAS'] +"ldsc_input/{model}_{trait}_{sp}_{h2}_{samplesize}.sumstats.gz",
          baselinef=expand(config_d['L2'] +"ukb_baseline_{{wind}}/{{model}}/{{source}}.{{samplesize}}.{chrom}.l2.ldscore.gz",chrom=range(1,23)),
          baseline_l2m=expand(config_d['L2'] +"ukb_baseline_{{wind}}/{{model}}/{{source}}.{{samplesize}}.{chrom}.l2.M",chrom=range(1,23)),
          baseline_l2m50=expand(config_d['L2'] +"ukb_baseline_{{wind}}/{{model}}/{{source}}.{{samplesize}}.{chrom}.l2.M_5_50",chrom=range(1,23))
      output:
          dataf="results/sim_ukb_{source}/{model}_{trait}_{sp}_{h2}_{samplesize}_{wind}.int.log"
      params:
          baseline=config_d['L2']+"ukb_baseline_{wind}/{model}/{source}.{samplesize}.",
          weights=config_d['L2']+"ukb_baseline_{wind}/{model}/{source}.{samplesize}.",#,config_d['WEIGHTS']+"weights.hm3_noMHC.",
          odir="results/sim_ukb_{source}/{model}_{trait}_{sp}_{h2}_{samplesize}_{wind}.int"
      shell:
          config_d['LDSC']+"ldsc.py --h2 {input.gwasf} --ref-ld-chr {params.baseline} --w-ld-chr {params.weights} --out {params.odir}"
```

`snakemake` sees that the rule `baseline_ldsc` can be used to create the file we asked for, because the output it generates: (`dataf="results/sim_ukb_{source}/{model}_{trait}_{sp}_{h2}_{samplesize}_{wind}.int.log") can be used to match `results/sim_ukb_ind/polym_3_10_0.5_10000_1.int.log` if it substitutes the following wildcards:

| wildcard     | substitution | meaning                                                                                                                                                        |   |
|--------------|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|---|
| {source}     | ind          | use in-sample LD scores                                                                                                                                        |   |
| {model}      | polym        | use the "polym" causal variant set.                                                                                                                            |   |
| {trait}      | 3            | Use the 3rd replicate of the gwas simulation for these parameters                                                                                              |   |
| {sp}         | 10           | Use the "10" subset of the "polym" causal variant set.  "10" represents 100% of the variants in "polym", "05" represents 50% of the variants and "01" is 12.5% |   |
| {h2}         | 0.5          | Use a gwas simulation with a total heritability of 0.5                                                                                                         |   |
| {samplesize} | 10000        | Use a gwas simulation with a sample size of 10000                                                                                                              |   |
| {wind}       | 1            | Use LD scores calculated with a 1 cM window.



After substituting wildcards, `snakemake` will then check to see if the following  files exist (note the `config_d` entries have been substituted out as well):

+ `"/gpfs/data/xhe-lab/polyg/" +"ldsc_input/polym_3_10_0.5_10000.sumstats.gz"`
+ `"/gpfs/data/xhe-lab/genomic_annotation/L2/" + "ukb_baseline_1/polym/ind.10000.{chrom}.l2.ldscore.gz"` (for  `chrom` in `1:22`)
+ `"/gpfs/data/xhe-lab/genomic_annotation/L2/" + "ukb_baseline_1/polym/ind.10000.{chrom}.l2.M"` (for `chrom` in `1:22`)
+ `"/gpfs/data/xhe-lab/genomic_annotation/L2/" + "ukb_baseline_1/polym/ind.10000.{chrom}.l2.M_5_50"` (for `chrom` in `1:22`)


If any of these files do not exist, the list of "rules" will be searched for a rule that has an `output:` section that matches the missing file (after substituting wildcards).  For example, if the file `/gpfs/data/xhe-lab/genomic_annotation/L2/ukb_baseline_1/polym/ind.10000.19.l2.M` was missing, it could be generated by the rule `cmp_baselineb_ldscores:` that looks like this:

file that matches the wil
```snakemake

  rule cmp_baselineb_ldscores:
      """ Compute baseline ld scores for the single annotation at a single chromosome"""
      input:
          famf=config_d["UKB_PED"]+"split_map/{source}_{samplesize}_chr{chrom}.fam",
          bimf=config_d["UKB_PED"]+"split_map/{source}_{samplesize}_chr{chrom}.bim",
          bedf=config_d["UKB_PED"]+"split_map/{source}_{samplesize}_chr{chrom}.bed",
          snp_list=config_d["MODELD"]+"model_{samplesize}/{model}.txt"
      output:
          l2=(config_d['L2']+"ukb_baseline_{wind}/{model}/{source}.{samplesize}.{chrom}.l2.M"),
          l2M_50=(config_d['L2']+"ukb_baseline_{wind}/{model}/{source}.{samplesize}.{chrom}.l2.M_5_50"),
          l2gz=(config_d['L2']+"ukb_baseline_{wind}/{model}/{source}.{samplesize}.{chrom}.l2.ldscore.gz")
      params:
          plink=config_d["UKB_PED"]+"split_map/{source}_{samplesize}_chr{chrom}",
          odir=config_d['L2']+"ukb_baseline_{wind}/{model}/{source}.{samplesize}.{chrom}",
          wind="{wind}"
      shell:
          config_d['LDSC']+"ldsc.py --l2 --bfile {params.plink}  --extract {input.snp_list} --ld-wind-cm {params.wind} --out {params.odir}"

```

If instead (or in addition), the file `/gpfs/data/xhe-lab/polyg/ldsc_input/polym_3_10_0.5_10000.sumstats.gz` were missing, the snakemake would see that the rule `sub_gwass_gcta_pm:` could create the file:

```snakemake
  rule sub_gwass_gcta_pm:
      input:
          vecf=config_d["GWAS"]+"gwas_covar_ss/{model}_{h2}_{trait}_{sp}_{samplesize}/res.fastGWA",
      output:
          tempf=config_d['GWAS'] +"ldsc_input/{model}_{trait}_{sp}_{h2}_{samplesize}.sumstats.gz"
      script:
          "../scripts/cols_cut_gcta.R"
```

The rule `sub_gwass_gcta_pm:`, rather than having a `shell:` section that gives the shell command for generating the output from the input, instead has a `script:` section. Looking in `PolygenicRSS/scripts/cols_cut_gcta.R` we see a very short R script:

```R
    library(dplyr)
    vroom::vroom(snakemake@input[["vecf"]]) %>% 
      dplyr::transmute(SNP=SNP,A1=A1,A2=A2,N=N,Z=BETA/SE) %>% 
    vroom::vroom_write(snakemake@output[["tempf"]],delim="\t")

```


In this R script the input file `vecf` is read in using the package `vroom` for quickly reading `tsv`s or `csv`s.  Z scores are calculated, then the data is written in `tsv` format to the file `tempf`. `snakemake` handles calling R in such a way that the `snakemake` object exists, and the `@input` and `@output` are filled correctly.


## Simulating GWAS

### Genotypes

To simulate GWAS you'll need genotypes, and if you're using the UK biobank you'll want to use a subset of the genotype data. For a gwas with 10000 individuals, you can run `snakemake   -j 400 --profile ../config/gardner "/scratch/${USER}/intersect_snplist/ukb_subset/all_10000_ind.bed"` ,where `${USER}` is your cri username, to generate the necessary plink-formatted genotype data that you'll need to run `simu`.

### Causal variants/genetic architecture

The simulation framework is currently only set up for simulating under a `GCTA`-style normal effect-size prior.  The file `"/scratch/${USER}/intersect_snplist/ukb_subset/model_{samplesize}/{model}.txt"` contains one variant-id per line, and every variant listed in the file will have a non-zero effect in the GWAS simulation. The file `"/scratch/${USER}/intersect_snplist/ukb_subset/model_{samplesize}/{model}.txt"`
