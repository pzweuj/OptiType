OptiType - Docker
========

Precision HLA typing from next-generation sequencing data

Authors: András Szolek, Benjamin Schubert, Christopher Mohr

Docker Update: pzweuj

Date: November 2024  

Version: 1.3.5.p

License: OptiType is released under a three-clause BSD license


Introduction
-------------
OptiType is a novel HLA genotyping algorithm based on integer linear
programming, capable of producing accurate 4-digit HLA genotyping predictions
from NGS data by simultaneously selecting all major and minor HLA Class I alleles.

[Original OptiType](https://github.com/FRED-2/OptiType)

Installation via Docker
-----------------------

1. Install Docker on your computer and make sure it works.

2. Call `docker pull ghcr.io/pzweuj/optitype:v1.3.5.p1 which will download the Docker image.

3. You can use the image as followes:

`docker run -v /path/to/data/folder:/data/ -t ghcr.io/pzweuj/optitype:v1.3.5.p1 -i input1 [input2] (-r|-d) -o /data/`

Usage
-------------

Optional step zero: you might want to filter your sequencing data for
HLA reads. Should you have to re-run OptiType multiple times on the same sample
(different settings, etc.) it could save you time if instead of giving OptiType
the full, multi-gigabyte sequencing data multiple times, you would rather give
it the relevant reads only, on the order of megabytes.

You can use any read mapper to do this step, although we suggest you use RazerS3.
Its only drawback is that due to way RazerS3 was designed, it loads all reads
into memory, which could be a problem on older, low-memory computing nodes.

Make sure to filter your files correctly depending on whether you have DNA
(exome, WGS) or RNA-Seq data. The reference fasta files are
`data/hla_reference_dna.fasta` and `data/hla_reference_rna.fasta` respectively.
Below is an example for DNA sequencing data:

```
>razers3 -i 95 -m 1 -dr 0 -o fished_1.bam /path/to/OptiType/data/hla_reference_dna.fasta sample_1.fastq

>samtools bam2fq fished_1.bam > sample_1_fished.fastq

>rm fished_1.bam
```

If you have paired-end data, repeat this with the second ends' fastq file as well.
Note: it's important that you filter the two ends individually. Don't use the
read mapper's paired-end capabilities.

After the optional filtering, OptiType can be called as follows:
```
>python /path/to/OptiTypePipeline.py -i sample_fished_1.fastq [sample_fished_2.fastq]
                    (--rna | --dna) [--beta BETA] [--enumerate N]
                    [-c CONFIG] [--verbose] --outdir /path/to/out_dir/
```

This will produce a time-stamped directory inside the specified output directory
containing a CSV file with the predicted optimal (and if enumerated, sub-optimal)
HLA genotype, and a pdf file containing a coverage plot of the predicted alleles
for diagnostic purposes.

```
>python OptiTypePipeline.py --help  

usage: OptiType [-h] --input FQ [FQ] (--rna | --dna) [--beta B]
                [--enumerate N] --outdir OUTDIR [--verbose] [--config CONFIG]

OptiType: 4-digit HLA typer

optional arguments:
  -h, --help            show this help message and exit
  --input FQ [FQ], -i FQ [FQ]
                        .fastq file(s) (fished or raw) or .bam files stored
                        for re-use, generated by an earlier OptiType run. One
                        file: single-end mode, two files: paired-end mode.
  --rna, -r             Use with RNA sequencing data.
  --dna, -d             Use with DNA sequencing data.
  --beta B, -b B        The beta value for for homozygosity detection (see
                        paper). Default: 0.009. Handle with care.
  --enumerate N, -e N   Number of enumerations. OptiType will output the
                        optimal solution and the top N-1 suboptimal solutions
                        in the results CSV. Default: 1
  --outdir OUTDIR, -o OUTDIR
                        Specifies the out directory to which all files should
                        be written.
  --prefix PREFIX, -p PREFIX
                        Specifies prefix of output files
  --verbose, -v         Set verbose mode on.
  --config CONFIG, -c CONFIG
                        Path to config file. Default: config.ini in the same
                        directory as this script
```

Furthermore, depending on your settings in `config.ini` you can choose to keep
the bam files OptiType produces when all-mapping reads against the reference:
these will be stored in the output directory of your current run.

Then, if you want to re-run OptiType on the same sample, you can give it those
intermediate `.bam` files as input instead of `.fastq` files, and spare on the
mapping part of the pipeline. Note: these `.bam` files have nothing to do with
those produced during the filtering/fishing step.


Test examples
-------------
DNA data (paired end):
```
python OptiTypePipeline.py -i ./test/exome/NA11995_SRR766010_1_fished.fastq ./test/exome/NA11995_SRR766010_2_fished.fastq --dna -v -o ./test/exome/
```

RNA data (paired end):
```
python OptiTypePipeline.py -i ./test/rna/CRC_81_N_1_fished.fastq ./test/rna/CRC_81_N_2_fished.fastq --rna -v -o ./test/rna/
```

Reference
-------------
Szolek, A, Schubert, B, Mohr, C, Sturm, M, Feldhahn, M, and Kohlbacher, O (2014).
OptiType: precision HLA typing from next-generation sequencing data
Bioinformatics, 30(23):3310-6.
