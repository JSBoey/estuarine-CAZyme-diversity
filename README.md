# CAZymes in an Estuary

There were some reshuffling of bins prior to cross-assembly dereplication (dRep). This updated analysis will use DASTool outputs for annotations.

## 0. Directory structure

Working directory: `/nesi/nobackup/uoa00348/boey/2024-estuarine-cazyme-diversity`

`data/` for inputs and outputs from automated processes.

`curation/` for any outputs that require manual curation. Outputs of manual curation need to have a space (column, section, etc) where the reason for selection is noted.

`db/` for newly obtained databases.

`bin/` for utility scripts and software.

`scripts/` and `slurm_out/` is self-evident.

## 1. Copy data and create query

Copy data

```sh
SRC=/nesi/nobackup/uoa00348/kim-waiwera/mags-unrefined

# Renamed DASTool bins
rsync --progress -av ${SRC}/00-mags-unrefined-rename-out data/
# DRAM outputs
rsync --progress -av ${SRC}/02-mags-unrefined-dram-out data/
# MAG scores
cp ${SRC}/03-read_mapping/summary_gtdbtk.drep.checkm.rename_unrefined.csv data/
```

Concatenate DRAM fasta files as HMMER and DIAMOND run better on a large file

```sh
(find data/ -name "genes.faa" -exec cat '{}' +) | \
    cut -f 1 -d ' ' \
    > data/dastool_bins.faa
```

Also create a concatenated genome and genes nucleotide fasta

```sh
find data/ -name "scaffolds.fna" -exec cat '{}' + > data/dastool_bins.scaffolds.fna
find data/ -name "genes.fna" -exec cat '{}' + > data/dastool_bins.genes.fna
```

Create a test query:

```sh
module load seqtk
seqtk sample -s1273 data/dastool_bins.faa 20000 > data/test_query.faa
```

## 2. Find CAZymes

**Download database**

```sh
# CAZy families
wget -P db/ https://bcb.unl.edu/dbCAN2/download/Databases/V13/dbCAN-HMMdb-V13.txt
mv db/dbCAN-HMMdb-V13.txt db/dbCAN-HMMdb-V13.hmm
wget -P db/ https://bcb.unl.edu/dbCAN2/download/Databases/fam-substrate-mapping-08012023.tsv
# Signal transduction
wget -P db/ https://bcb.unl.edu/dbCAN2/download/Databases/V12/stp.hmm
# Transcription factors
wget -P db/ https://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-1.hmm
wget -P db/ https://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-2.hmm
# CAZyDB
wget -P db/ https://bcb.unl.edu/dbCAN2/download/Databases/V13/CAZyDB.07142024.fa
mv db/CAZyDB.07142024.fa db/CAZyDB.07142024.faa
# Transporters
wget -P db/ http://www.tcdb.org/public/tcdb
mv db/tcdb db/tcdb_20241918.faa
# CAZyme subfamilies
wget -P db/ https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN_sub.hmm
mv db/dbCAN_sub.hmm db/dbCAN_sub_20220830.hmm
wget -P db/ https://bcb.unl.edu/dbCAN2/download/Databases/subfam_EC_mapping.tsv
mv db/subfam_EC_mapping.tsv db/subfam_EC_mapping_20240105.tsv
wget -P db/ https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-sub.substrate.mapping.xls
# PUL/CGC
wget -P db/ https://bcb.unl.edu/dbCAN2/download/Databases/PUL.faa
mv db/PUL.faa db/PUL_20240105.faa
wget -P db/ https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_12-12-2023.xlsx
```

**Annotate against databases**

`scripts/hmmer_dbs.sl`

```sh
#!/bin/bash -e
#SBATCH --job-name=annot_hmmer
#SBATCH --account=ga02676
#SBATCH --time=24:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err
#SBATCH --array=0-4
#SBATCH --partition=milan

# Modules
module purge
module load HMMER/3.3.2-GCC-12.3.0

# Array variables
a=(db/*.hmm)
db=${a[$SLURM_ARRAY_TASK_ID]}
db_base=$(basename ${db} .hmm)

seqin=data/dastool_bins.faa
domtblout=data/annotations/dastool_bins.${db_base}.domtblout
out=${domtblout/domtbl/hmm}
z=$(hmmstat ${db} | grep -cv '^#')

if [[ ! -d $(dirname ${domtblout}) ]]; then
    mkdir -p $(dirname ${domtblout})
fi

echo "
This is ${SLURM_JOB_ID} searching ${seqin} against ${db_base} (Z = ${z})
"

hmmsearch --cpu $SLURM_CPUS_PER_TASK \
          --domtblout $domtblout \
          -o $out \
          -Z $z \
          $db \
          $seqin
```

JobID: 50271991 (dbCAN-sub didn't finish)
JobID: 50281972

`scripts/annot_dmnd.sl`

```sh
#!/bin/bash -e
#SBATCH --job-name=annot_dmnd
#SBATCH --account=ga02676
#SBATCH --time=24:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err
#SBATCH --array=0-2
#SBATCH --partition=milan

# Modules
module purge
module load DIAMOND/2.1.9-GCC-11.3.0

# Array variables
a=(db/*.faa)
db=${a[$SLURM_ARRAY_TASK_ID]}
db_base=$(basename ${db} .faa)

seqin=data/dastool_bins.faa
out=data/annotations/dastool_bins.${db_base}.b6o
evalue=1e-10
outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qcovhsp scovhsp bitscore corrected_bitscore evalue'

if [[ ! -d $(dirname ${out}) ]]; then
    mkdir -p $(dirname ${out})
fi

# Make database if doesn't exist
if [[ ! -f ${db/faa/dmnd} ]]; then
    echo "Need to build database..."
    diamond makedb --threads $SLURM_CPUS_PER_TASK --db ${db/faa/dmnd} --in ${db}
fi

# CAZymes need more stringent E-value
if echo ${db_base} | grep -vq "tcdb"; then
    evalue=1e-102
fi

diamond blastp --threads $SLURM_CPUS_PER_TASK --db ${db/faa/dmnd} \
    --evalue ${evalue} --query ${seqin} \
    --out ${out} --outfmt ${outfmt}
```

JobID: 50271982 (Complete)

Forgot to add KEGG... Ugh...

Format KEGG database

```bash
module purge
module load DIAMOND/2.1.9-GCC-11.3.0

dbdir=/nesi/project/uoa02469/Databases/KEGG_2024

diamond makedb --in ${dbdir}/refprok.pep.gz --db db/kegg_2024.refprok.dmnd --threads 12
zcat ${dbdir}/refprok.dat.gz > db/kegg_2024.refprok.dat
```

`scripts/dmnd_kegg.sl`

```sh
#!/bin/bash -e
#SBATCH --job-name=dmnd_kegg
#SBATCH --account=ga02676
#SBATCH --time=3:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err
#SBATCH --partition=milan

# Modules
module purge
module load DIAMOND/2.1.9-GCC-11.3.0

db=db/kegg_2024.refprok.dmnd
db_base=$(basename ${db} .dmnd)

seqin=data/dastool_bins.faa
out=data/annotations/dastool_bins.${db_base}.b6o

evalue=1e-102
outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qcovhsp scovhsp bitscore corrected_bitscore evalue'

if [[ ! -d $(dirname ${out}) ]]; then
    mkdir -p $(dirname ${out})
fi

diamond blastp --threads $SLURM_CPUS_PER_TASK --db ${db} \
    --evalue ${evalue} --query ${seqin} \
    --out ${out} --outfmt ${outfmt}
```

JobID: 50385522

## 3. Find domains

Using InterProScan 5.69-101.0 (released 25 July 2024)

Download and unpack

```sh
cd bin/

wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.69-101.0/interproscan-5.69-101.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.69-101.0/interproscan-5.69-101.0-64-bit.tar.gz.md5

md5sum -c interproscan-5.69-101.0-64-bit.tar.gz.md5

tar -pxvzf interproscan-5.69-101.0-*-bit.tar.gz
```

Setup and test

```sh
cd interproscan-5.69-101.0

python3 setup.py -f interproscan.properties

./interproscan.sh -i test_all_appl.fasta -f xml

ln -sr ./interproscan.sh ../
```

Preformat and split inputs (InterProScan is slow...)

```sh
mkdir -p ./tmp/iprin

module load SeqKit/2.4.0
cat data/dastool_bins.faa | \
    sed -e 's/\*//g' | \
    seqkit split2 --by-size 80000 --out-dir tmp/iprin \
                  --by-size-prefix dastool_bins.part_
```

`scripts/ipr5.sl`

```sh
#!/bin/bash -e
#SBATCH --job-name=ipr5
#SBATCH --account=ga02676
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-35
#SBATCH --partition=milan

# Modules
module load Java/11.0.4

# Add to path
export PATH=$PATH:$(pwd)/bin

# Array variables
a=(tmp/iprin/*)
input=${a[$SLURM_ARRAY_TASK_ID]}
inbase=$(basename ${input} .fasta)
outbase=tmp/iprout/${inbase}.ipr5
tempdir=./tmp/iprtmp/${SLURM_JOB_ID}

# Run variables
outfmt=xml,tsv
exclappl='ProSiteProfiles,ProSitePatterns' # Need PCRE2 library
cpu=$(( SLURM_CPUS_PER_TASK - 2 ))

if [[ ! -d ${tempdir} ]]; then
    mkdir -p ${tempdir}
fi

if [[ ! -d $(dirname ${outbase}) ]]; then
    mkdir -p $(dirname ${outbase})
fi

interproscan.sh \
    --input ${input} \
    --output-file-base ${outbase} \
    --formats ${outfmt} \
    --excl-applications ${exclappl} \
    --cpu ${cpu} \
    --tempdir ${tempdir} \
    --pathways --goterms
```

JobID: 50351954

There are random failures... Check which ones failed:

```bash
sacct -j50351954 --state "FAILED" --format jobid | \
    tail -n+3 | \
    grep -v "[\.\+]"
```

Arrays 9, 14, 16, 19-21, 25, 28, 30 didn't work...

Change array ID in script and re-run

JobID: 50381587

----

I'm not entirely happy with the results obtained from a BLASTp of TCDB. During CGC finding, I found some genes that were more likely than not to be part of a CGC, but have very distant homology to known transporters (excellent E-value and bitscores, < 30% ID). I'll try a different tactic: Include sequences that match TCDoms.

```sh
cd db/

wget https://tcdb.org/public/tcDoms.tar.gz

tar -xzvf tcDoms.tar.gz
```

`scripts/annot_tcDoms.sl`

```bash
#!/bin/bash -e
#SBATCH --job-name=tcdoms
#SBATCH --account=ga02676
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%j.%x.out
#SBATCH --error=slurm_err/%j.%x.err
#SBATCH --partition=milan

# Modules
module purge
module load HMMER/3.3.2-GCC-12.3.0

# Variables
db=db/tcDoms/tcDomsGlobal/tcDoms.hmm
seqin=data/dastool_bins.faa
domtblout=data/annotations/dastool_bins.tcdoms.domtblout
out=${domtblout/domtbl/hmm}
z=$(hmmstat ${db} | grep -cv '^#')

if [[ ! -d $(dirname ${domtblout}) ]]; then
    mkdir -p $(dirname ${domtblout})
fi

echo "
This is ${SLURM_JOB_ID} searching ${seqin} against tcDoms.hmm (Z = ${z})
"

# Run
hmmsearch --cpu $SLURM_CPUS_PER_TASK \
          --domtblout $domtblout \
          -o $out \
          -Z $z \
          $db \
          $seqin

```

JobID: 50568035

## 3a. Find dbCAN-sub domains

The HMMER3 script using the entire query timed out with increased time to 24 hours. I'll use the InterProScan input (80K seqs per file) for dbCAN-sub annotations in an array.

`scripts/annot_hmmer.dbcansub.array.sl`

```sh
#!/bin/bash -e
#SBATCH --job-name=dbcansub
#SBATCH --account=ga02676
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%j.%A_%a.%x.out
#SBATCH --error=slurm_err/%j.%A_%a.%x.err
#SBATCH --array=0-35
#SBATCH --partition=milan

# Modules
module purge
module load HMMER/3.3.2-GCC-12.3.0

# Array variables
a=(tmp/iprin/*)
seqin=${a[$SLURM_ARRAY_TASK_ID]}
inbase=$(basename ${seqin} .fasta)

db=db/dbCAN_sub.hmm
domtblout=tmp/dbcansubout/${inbase}.domtblout
out=${domtblout/domtbl/hmm}
z=$(hmmstat ${db} | grep -cv '^#')

if [[ ! -d $(dirname ${domtblout}) ]]; then
    mkdir -p $(dirname ${domtblout})
fi

echo "
This is ${SLURM_JOB_ID} searching ${seqin} against ${db_base} (Z = ${z})
"

hmmsearch --cpu $SLURM_CPUS_PER_TASK \
          --domtblout $domtblout \
          -o $out \
          -Z $z \
          $db \
          $seqin
```

JobID: 50351297

## 4. Parse outputs

For CAZy predictions by HMMER:

`bin/hmmsearch_parser.sh`

This is a modification of the original parser script `hmmscan_parser.sh` by Yanbin Yin.

```sh
#!/bin/bash -e

# Modified from hmmscan-parser.sh
# Originally written by Yanbin Yin
# Usage: hmmsearch-parser.sh hmmsearch_domain_table
# Removed filtering for E-value and domain coverage to retain flexibility

cat $1 | \
    grep -v '^#' | \
    awk '{ print $4, $6, $1, $3, $13, $16, $17, $18, $19 }' | \
    sed 's/ /\t/g' | \
    sort -k 3,3 -k 8n -k 9n | \
    perl -e '
    while (<>) {
        chomp;
        @a = split;
        next if($a[-1] == $a[-2]);
        push(@{ $b{ $a[2] } }, $_);
    }

    foreach (sort keys %b) {
        @a = @{ $b{$_} };
        for ($i = 0; $i < $#a; $i++) {
            @b = split(/\t/, $a[$i]);
            @c = split(/\t/, $a[ $i + 1 ]);
            $len1 = $b[-1] - $b[-2]; # Aligned length b
            $len2 = $c[-1] - $c[-2]; # Aligned length c
            $len3 = $b[-1] - $c[-2]; # Overlap between b and c
            if($len3 > 0 
                and ($len3 / $len1 > 0.5 or $len3 / $len2 > 0.5)) 
            {
                if ($b[4] < $c[4]) { # Remove lower E-value if overlap
                    splice(@a, $i + 1, 1);
                } else {
                    splice(@a, $i, 1);
                }
                $i = $i - 1;
            }
        }
        foreach(@a) { print $_."\n"; }
    }
    ' | \
    perl -e '
    while(<>) {
        chomp;
        @a = split(/\t/, $_);
        if (($a[-1] - $a[-2]) > 80) {
            print $_, "\t", ($a[-3] - $a[-4]) / $a[1], "\n" if $a[4] < 1e-5;
        } else {
            print $_, "\t", ($a[-3] - $a[-4]) / $a[1], "\n" if $a[4] < 1e-3;
        }
    }
    ' | \
    awk '$NF > 0.30' | \
    sort -k 3 -k 8,9g | \
    awk '
    BEGIN {
        OFS="\t"
    }
    { print $3, $4, $1, $2, $5, $6, $7, $8, $9, $10 }
    ' | \
    sed '1 i\query\tquery_len\ttarget\ttarget_length\tevalue\ttarget_from\ttarget_to\tquery_from\tquery_to\tdomain_coverage'
```

For non-dbCAN_sub CAZyme predictions:

```sh
cd data/annotations

for i in *.domtblout; do
    echo "Parsing ${i}"
    hmmsearch_parser.sh ${i} > ${i}_parsed
done

cd ../../
```

For dbCAN_sub predictions:

The concatenated `.domtblout` files are huge. Parse individual files then merge parsed outputs.

`script/parse_dbcansub.sl`

```sh
#!/bin/bash -e
#SBATCH --job-name=parse_dbcansub
#SBATCH --account=ga02676
#SBATCH --time=10:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm_out/%j.%A_%a.%x.out
#SBATCH --error=slurm_err/%j.%A_%a.%x.err
#SBATCH --array=0-35
#SBATCH --partition=milan

export PATH=$PATH:$(pwd)/bin

a=(tmp/dbcansubout/*.domtblout)
b=${a[$SLURM_ARRAY_TASK_ID]}
out=${b}_parsed

echo "Parsing ${b}"

hmmsearch_parser.sh ${b} > ${out}
```

JobID: 50359268

When finished

```sh
a=(tmp/dbcansubout/*_parsed)
pout=data/annotations/dastool_bins.dbCAN_sub.domtblout_parsed

head -n1 ${a[0]} > ${pout}
for i in ${a[@]}; do
    echo "$(basename ${i}) >> $(basename ${pout})"
    tail -n+2 ${i} >> ${pout}
done

# Check
b=0
for i in ${a[@]}; do
    n=$(tail -n+2 ${i} | wc -l | cut -f1 -d' ')
    b=$((b + n))
done
echo $((b+1))
```

Archive and compress partitioned outputs

```sh
chmod 777 tmp/dbcansubout/*
tar -cvf data/annotations/dastool_bins.dbCAN_sub.all_outputs.tar tmp/dbcansubout/*

module load pigz
pigz -p 8 data/annotations/dastool_bins.dbCAN_sub.all_outputs.tar
```

For diamond outputs:

```sh
cd data/annotations

ofmt="qseqid sseqid pident length mismatch gapopen qstart qend sstart send qcovhsp scovhsp bitscore corrected_bitscore evalue"
hfmt=$(echo ${ofmt} | sed 's/ /\t/g')

for i in *.b6o; do
    sed "1i ${hfmt}" ${i} > ${i}_parsed
done
```

## 4a. Parsing InterProScan

Archive header-filled parts to `data/annotations`

```bash
cd tmp/iprout

module load pigz

ofmt="query query_md5 query_length analysis signature_accession signature_description query_start query_stop signature_score match_status run_date interpro_accession interpro_description go_terms pathways"
hfmt=$(echo ${ofmt} | sed 's/ /\t/g')

for i in *.tsv; do
  sed "1i ${hfmt}" ${i} > ${i/ipr5/ipr5_parsed}
done

tar -cvf ../../data/annotations/dastool_bins.ipr5_parsed_tsv.tar *_parsed.tsv
tar -cvf ../../data/annotations/dastool_bins.ipr5_xml.tar *.xml

cd ../../data/annotations

pigz -vp8 dastool_bins.ipr5_parsed_tsv.tar
pigz -vp8 dastool_bins.ipr5_xml.tar
tar -cvf data/annotations/dastool_bins.ipr5.all_tsv.tar tmp/iprout/*.tsv
pigz -vp8 data/annotations/dastool_bins.ipr5.all_tsv.tar
```

The size of the data is difficult to move. Need to convert the TSV files into Parquet.

```bash
module load R/4.3.2-foss-2023a
R
```

```R
library(arrow)
library(dplyr)

datafiles <- list.files("iprout", pattern=".*.tsv", full.names=T)
dts <- open_dataset(datafiles, format="tsv")

set_cpu_count(num_threads=12)
write_dataset(dts, "ipr5_parsed_parquet", 
              format="parquet", 
              partitioning="analysis")
```

## 5. Matching metadata

Based on CAZy's webpages

```r
#!/usr/bin/env Rscript

# Downloads metadata for characterized CAZymes from www.cazy.org.
# Example usage:
#   Rscript download_cazy_tables.R -i cazy_families.txt

# Packages  ----
pkgs <- c("XML", "tidyverse", "furrr", "progressr", "optparse")
suppressPackageStartupMessages({
  for (p in pkgs) {
    if (!require(p, character.only = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org/")
    }
  
    library(p, character.only = TRUE)
  }
})

# Functions ----
# Create HTML link based on CAZyme family input
makeAddress <- function(cazyme_family) {
  paste0("http://www.cazy.org/", cazyme_family, "_characterized.html")
}

# Are tables paginated?
# If yes, function will return vector of page links
# Otherwise, return an empty vector
makePageLinks <- function(html_link) {
  # Read HTML as text lines
  html_vector <- readLines(html_link)
  
  # Find pagination lines
  match_vector <- str_subset(html_vector, '\\?debut_FUNC')
  
  # If there are HTML lines indicating pagination, create list of links, otherwise return
  # HTML link
  if (is_empty(match_vector)) {
    result <- html_link
  } else {
    page_link <- unique(
      unlist(
        str_extract_all(match_vector, "\\?debut_FUNC=\\d+#pagination_FUNC")
      )
    )
    page_link_vector <- paste0(html_link, page_link)
    result <- append(html_link, page_link_vector)
  }
  
  return(result)
}

# Download and concatenate tables
downloadTables <- function(page_links) {
  # Download HTML and parse into table
  tbs <- map(page_links, ~ {
    html_lines <- readLines(.x)
    
    # Replace breaks '<br>' with spaces
    parsed_breaks <- str_replace_all(html_lines, "<br>", " ")
    
    # Parse downloaded HTML lines
    parsed_html <- readHTMLTable(parsed_breaks)[[2]]
  })
  
  # Bind rows
  if (length(tbs) > 1) {
    bind_rows(tbs)
  } else {
    tbs[[1]]
  }
}

# Re-format tables
formatTable <- function(tb) {
  headers <- c("protein_name", "ec_number", "reference", "organism", "genbank", "uniprot", "pdb3d")
  
  tb0 <- tb %>% 
    # Remove unnecessary rows
    filter(!str_detect(V1, "suivante|Protein Name|Top")) %>% 
    # Clean up cells
    mutate(
      across(everything(), ~ str_trim(.x))
    )
  
  # Append headers
  if (ncol(tb0) > 7) {
    colnames(tb0) <- append(headers, "subfamily")
  } else {
    colnames(tb0) <- headers
  }
  
  return(tb0)
}

# Append taxonomic domain as column
appendTaxaDomain <- function(formatted_table) {
  # Create vector of where to split
  split_indicator <- is.na(formatted_table[["organism"]])
  
  # Split table by taxonomic domain
  split_table <- split(formatted_table, cumsum(split_indicator))
  
  # Get taxonomic domains
  taxa_domain <- map_chr(split_table, ~ .x[1, 1])
  split_table <- set_names(split_table, taxa_domain)
  taxa <- unique(taxa_domain) %>% 
    set_names(., .)
  
  # Append taxonomic domain to each row
  domain_tables <- map(taxa, ~ {
    keep(split_table, names(split_table) == .x) %>% 
      bind_rows() %>% 
      filter(!is.na(organism))
  }) %>% 
    bind_rows(.id = "taxa_domain")
}

# Wrapper function ----
getCAZyTables <- function(cazyme_family) {
  
  if (!is.character(cazyme_family) || length(cazyme_family) > 1) {
    stop("Input must be a character vector of length 1.")
  }
  
  output <- makeAddress(cazyme_family) %>%
    makePageLinks() %>%
    downloadTables() %>%
    formatTable() %>%
    appendTaxaDomain()
  
  gc()
  
  output
}

getTablesSafely <- possibly(getCAZyTables, otherwise = "Not found")

# Parse options ----
optlist <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "File with CAZyme families (one entry per line; '-' for stdin)"),
  make_option(c("-o", "--outfile"), type = "character", 
              default = "./characterized_cazy.tsv",
              help = "Output file"),
  make_option(c("-t", "--threads"), type = "integer", 
              default = 1,
              help = "Number of parallel processes [default = 1]")
)

prog_description <- "
Download and parse tables of characterized CAZymes from CAZy database.
"

opt <- parse_args(
  OptionParser(option_list = optlist, 
               description = prog_description)
)

# Validate ----
if (is.null(opt$input))
  stop("Input file required")

if (opt$threads < 1)
  stop("-t must be a positive integer")

# Run ----
if (opt$input == "-") {
  cfam <- readLines(file("stdin"))
} else {
  cfam <- readLines(opt$input)
}

cfam <- str_subset(cfam, "AA|CBM|CE|GH|GT|PL")
names(cfam) <- cfam

cat(str_glue("Found {length(cfam)} CAZyme families"), "\n")

cat(str_glue("Downloading tables with {opt$threads} threads"), "\n")

if (opt$threads > 1L) {
  plan(multisession, workers = opt$threads)
}

with_progress({
  p <- progressor(steps = length(cfam))
  clist <- future_map(cfam, function(x) {
    p()
    getTablesSafely(x)
  })
})

plan(sequential)

# Clean up table
cat("Cleaning table\n")

ctab <- keep(clist, is.data.frame) %>%
  bind_rows(.id = "family")

write_tsv(ctab, opt$outfile)

cat(str_glue("Wrote {opt$outfile} containing {nrow(ctab)} entries"), "\n\n")
```

```sh
cd db/

module load R/4.3.2-foss-2023a

R_PROGRESSR_ENABLE=TRUE
grep '>' CAZyDB.07142024.faa | \
    sed -r 's/.*\|([A-Z]+[0-9]+).*/\1/g' | \
    sort -u | \
    download_cazy_tables.R -i - -t 10
```

## 6. Representative genomes

These genomes are dereplicated at 98% ANI with pre-calculated CheckM values.

```bash
cp /nesi/nobackup/uoa00348/kim-waiwera/mags-unrefined/01-mags-unrefined-checkm-v121-out/genomes_checkm.txt data/
```

CheckM output needs to be in a specific format:

```text
genome,completeness,contamination,strain_heterogeneity
<genome1>,<completeness1>,<contamination1>,<strain_heterogeneity1>
<genome2>,<completeness2>,<contamination2>,<strain_heterogeneity2>
<genome3>,<completeness3>,<contamination3>,<strain_heterogeneity3>
```

```bash
cd data
echo "genome,completeness,contamination,strain_heterogeneity" > checkm_for_drep.csv
tail -n+2 genomes_checkm.txt | \
  cut -f1,12-14 | \
  awk -v OFS=, '{ print $1".fa", $2, $3, $4 }' \
  >> checkm_for_drep.csv
cd ../
```

Run dRep

`scripts/drep_98ani.sl`

```bash
#!/bin/bash -e
#SBATCH --job-name=drep_98ani
#SBATCH --account=ga02676
#SBATCH --time=3:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=24
#SBATCH --output=slurm_out/%j.%x.out
#SBATCH --error=slurm_err/%j.%x.err
#SBATCH --partition=milan

module purge
module load drep/3.4.2-gimkl-2022a-Python-3.10.5

genome_dir=data/00-mags-unrefined-rename-out
work_dir=data/drep_98ani_out
checkm_csv=data/checkm_for_drep.csv

dRep dereplicate \
  --processors ${SLURM_CPUS_PER_TASK} \
  --genomes ${genome_dir}/*.fa \
  --genomeInfo ${checkm_csv} \
  --S_algorithm gANI \
  --completeness 0 \
  --contamination 100 \
  --P_ani 0.90 \
  --S_ani 0.98
```

Job ID: 50683528

Compile list of winning genomes that also have >= 50% completeness and <= 5% contamination after considering strain heterogeneity.

```bash
cd data/
mkdir -p drep_98ani_50comp_5cont

ls -1 drep_98ani_out/dereplicated_genomes/*.fa | \
  sed 's/.*\///g' | \
  # Winning dRep 98% ANI genomes
  grep -f - checkm_for_drep.csv | \
  # Contamination after accounting strain heterogenety
  awk -v FS=, -v OFS=, '$2 >= 50 && $3*(100-$4)*.01 <= 5' | \
  tail -n+2 | \
  cut -f1 -d, | \
  xargs -I % cp --verbose drep_98ani_out/dereplicated_genomes/% drep_98ani_50comp_5cont/

# Append bin name to FASTA headers
cd drep_98ani_50comp_5cont

for i in *.fa; do
  base=$(basename ${i} .fa)
  echo "Appending to ${base}"
  sed -i "s/>/>${base}_/g" ${i}
done

cd ../../
```

## 7. Differential coverage

Copy read libraries bins into `data/`

WGS data is in `/nesi/project/ga02676/Waiwera_project/1.Qual_filtered_trimmomatic/`
WTS data is in `/nesi/project/ga02676/Waiwera_project/sze_works/1.WTS/3.trimmomatic/`

```bash
cd data/

mkdir w{g,t}s
rsync -av --progress /nesi/project/ga02676/Waiwera_project/1.Qual_filtered_trimmomatic/*_R{1,2}* wgs/
rsync -av --progress /nesi/project/ga02676/Waiwera_project/sze_works/1.WTS/3.trimmomatic/*.fastq.gz wts/
```

Water samples for WGS were run on separate lanes, concatenate before proceeding.

```bash
cd wgs

for i in {1..9}; do
  for j in {1..2}; do
    echo "Appending Filt.S${i}_L?_R${j}.fastq.gz"
    cat Filt.S${i}_L?_R${j}.fastq.gz > Filt.S${i}_R${j}.fastq.gz
  done
done

# Remove individual lane files
rm Filt.S?_L?_R?.fastq.gz
```

Files have different naming conventions, rename to the following:

`<read_type>.<sample_type>.<site>_<biological_replicate>.<read_pair>.fastq.gz`

```bash
# Still in data/wgs/
for i in Filt*; do
  read_type=WGS
  sample_type=Filt
  site=$(echo $i | sed -E 's/.*S([0-9]).*/\1/g')
  biol_rep=1
  read_pair=$(echo $i | sed -E 's/.*(R[0-9]).*/\1/g')
  newname=${read_type}.${sample_type}.S${site}_${biol_rep}.${read_pair}.fastq.gz
  mv --verbose $i $newname
done

for i in Sed*; do
  read_type=WGS
  sample_type=Sed
  site=$(echo $i | sed -E 's/.*S([0-9]).*/\1/g')
  biol_rep=$(echo $i | sed -E 's/.*Sample([0-9]).*/\1/g')
  read_pair=$(echo $i | sed -E 's/.*(R[0-9]).*/\1/g')
  newname=${read_type}.${sample_type}.S${site}_${biol_rep}.${read_pair}.fastq.gz
  mv --verbose $i $newname
done

cd ../wts

for i in *Filt*; do
  read_type=WTS
  sample_type=Filt
  site=$(echo $i | sed -E 's/.*Filt([0-9]).*/\1/g')
  biol_rep=1
  read_pair=$(echo $i | sed -E 's/.*(R[0-9]).*/\1/g')
  newname=${read_type}.${sample_type}.S${site}_${biol_rep}.${read_pair}.fastq.gz
  mv --verbose $i $newname
done

for i in *Sed*; do
  read_type=WTS
  sample_type=Sed
  site=$(echo $i | sed -E 's/.*S([0-9]).*/\1/g')
  biol_rep=$(echo $i | sed -E 's/.*Sed([0-9]).*/\1/g')
  read_pair=$(echo $i | sed -E 's/.*(R[0-9]).*/\1/g')
  newname=${read_type}.${sample_type}.S${site}_${biol_rep}.${read_pair}.fastq.gz
  mv --verbose $i $newname
done
```

Index genomes

```bash
cd data/

mkdir -p read_mapping/{index,bam}

module purge
module load Bowtie2/2.5.4-GCC-12.3.0

bowtie2-build --threads 12 \
  drep_98ani_50comp_5cont/*.fa \
  read_mapping/index/drep_98ani_50comp_5cont
```

Figure out rough insert size across samples based on a subset of reads

`scripts/insert_size.sl

```bash
#!/bin/bash -e
#SBATCH --job-name=insize
#SBATCH --account=uoa00348
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-65
#SBATCH --partition=milan

a=(data/w{g,t}s/*.R1.fastq.gz)
r1=${a[$SLURM_ARRAY_TASK_ID]}
r2=${r1/R1/R2}
base=$(basename ${r1} .R1.fastq.gz)
outdir=data/read_mapping/insert_size_estimates

# Subsample reads
module purge
module load seqtk/1.4-GCC-11.3.0

sub_n=1000000
seed=1647

seqtk sample -s ${seed} ${r1} ${sub_n} > ${outdir}/${base}.R1.fq
seqtk sample -s ${seed} ${r2} ${sub_n} > ${outdir}/${base}.R2.fq

# Map reads
module purge
module load Bowtie2/2.5.4-GCC-12.3.0 SAMtools/1.19-GCC-12.3.0

index=$(dirname ${outdir})/index/drep_98ani_50comp_5cont

bowtie2 --sensitive -x ${index} -p $SLURM_CPUS_PER_TASK -1 ${r1} -2 ${r2} | \
  samtools sort -@ $SLURM_CPUS_PER_TASK -o ${outdir}/${base}.bam

# Get metrics
samtools stats -@ $SLURM_CPUS_PER_TASK ${outdir}/${base}.bam > ${outdir}/${base}.samstats

module purge
module load R/4.3.2-foss-2023a picard/2.26.10-Java-11.0.4

picard CollectInsertSizeMetrics \
  I=${outdir}/${base}.bam \
  O=${outdir}/${base}.insert_size_metrics.txt \
  H=${outdir}/${base}.insert_size_histogram.pdf

```

Job ID: 50694521

PICARD insert size estimates:

```bash
sed -n '7p' WGS.Filt.S1_1.insert_size_metrics.txt | \
    sed 's/^/SAMPLE\t/' > picard.insert_size_metrics.txt
for i in W{G,T}S.*.insert_size_metrics.txt; do
    base=$(basename ${i} .insert_size_metrics.txt)
    echo "Appending ${base}"
    metric=$(sed -n '8p' ${i})
    echo -e "${base},${metric}" >> picard.insert_size_metrics.txt
done
sed -i 's/\t/,/g' picard.insert_size_metrics.txt
```

A reasonable estimate for `-I` and `-X`:

|     | `-I` | `-X` |
| --- | ---- | ---- |
| WGS | 200  | 800  |
| WTS |   0  | 500  |

The setting for WTS is default.

Map reads

`script/bt2map.sl`

```bash
#!/bin/bash -e
#SBATCH --job-name=bt2map
#SBATCH --account=uoa00348
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-65
#SBATCH --partition=milan

# Read mapping
module purge
module load Bowtie2/2.5.4-GCC-12.3.0 SAMtools/1.19-GCC-12.3.0

a=(data/w{g,t}s/*.R1.fastq.gz)
r1=${a[$SLURM_ARRAY_TASK_ID]}
r2=${r1/R1/R2}
base=$(basename ${r1} .R1.fastq.gz)
outdir=data/read_mapping/bam
index=$(dirname ${outdir})/index/drep_98ani_50comp_5cont
minin=0
maxin=500

# Modify insert size parameters
if grep -q "WGS" ${base}; then
    minin=200
    maxin=800
fi

bowtie2 --sensitive \
    -x ${index} \
    -p ${SLURM_CPUS_PER_TASK} \
    -I ${minin} -X ${maxin} \
    -1 ${r1} -2 ${r2} | \
    samtools sort -@ $SLURM_CPUS_PER_TASK \
        -o ${outdir}/${base}.bam -
```

JobID: 50748118

Get mapping metrics

```bash
cd data/read_mapping/bam

module purge
module load SAMtools/1.19-GCC-12.3.0

for i in *.bam; do
    echo ${i}
    samtools flagstat -@8 -Otsv ${i} > ${i/.bam/.stat}
done

for i in *.stat; do
    base=$(basename ${i} .stat)
    echo "Appending ${base}"
    sed "s/^/${base}\t/g" ${i} >> ../alignment_statistics.tsv
done

cd ../../
```

Get contig coverage information

```bash
module purge
module load CoverM/0.7.0-GCC-12.3.0

bamfiles=$(ls bam/*.bam | tr "\n" " ")
coverm contig \
    --bam-files ${bamfiles} \
    --methods length covered_bases count \
    --threads 8 \
    --output-file cov/contig.len_base_count.tsv

coverm contig \
    --bam-files ${bamfiles} \
    --methods metabat \
    --threads 8 \
    --output-file cov/contig.metabat_cov.tsv

cd ../
```

Gather gene features for feature counting

```bash
saf=read_mapping/cov/drep_98ani_50comp_5cont.saf
printf "%s\t%s\t%s\t%s\t%s\n" GeneID Chr Start End Strand > ${saf}
for i in drep_98ani_50comp_5cont/*.fa; do
    base=$(basename ${i} .fa)
    echo "Gathering features from ${base}/genes.gff"
    sed 's/;/\t/g' 02-mags-unrefined-dram-out/${base}/genes.gff \
        | grep -v '^#' \
        | cut -f1,4,5,7,9 \
        | sed 's/ID=//g' \
        | awk -v FS="\t" -v OFS="\t" '{print $5, $1, $2, $3, $4}' \
        >> ${saf}
done
```

Count features

```bash
cd data/read_mapping/bam

module purge
module load Subread/2.0.7-GCC-12.3.0

# Count fragments
featureCounts -p -T 12 -F SAF \
    --countReadPairs \
    -a ../cov/drep_98ani_50comp_5cont.saf \
    -o ../cov/WTS_featureCounts.fragments.tsv \
    WTS*.bam
```

## 8. PULs

Okay... So... This late in the game I find out dbCAN3 uses BLASTp to get PULs to match with lower thresholds that expected... I get it, BLASTp is waay more sensitive, but I shouldn't have to dig through code to find it... 

I'm not going to bother subsetting the data to only CGCs... that might change but at least if I blast against the whole set I won't be missing things

Set up database

```bash
cd db/

module purge && module load BLAST/2.16.0-GCC-12.3.0

makeblastdb -dbtype prot \
            -in PUL_20240105.faa \
            -input_type fasta \
            -out PUL_20240105_blastdb
```

`scripts/pul_blastp.sl`

```bash
#!/bin/bash -e
#SBATCH --job-name=pul_blastp
#SBATCH --account=uoa00348
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%x.%j.%x.out
#SBATCH --error=slurm_err/%x.%j.%x.err
#SBATCH --partition=milan

module purge && module load BLAST/2.16.0-GCC-12.3.0

db=db/PUL_20240105_blastdb
query=data/dastool_bins.faa
evalue=0.01
outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp"
out=data/annotations/dastool_bins.PUL_20240105.blastp.b6o

blastp -task blastp -query ${query} -max_hsps 1 -db ${db} -outfmt "${outfmt}" -out ${out} -evalue ${evalue} -num_threads ${SLURM_CPUS_PER_TASK}

header=$(echo ${outfmt} | sed -e 's/6 //g' -e 's/ /\t/g')
sed -i "1i${header}" ${out}

awk -v FS="\t" -v OFS="\t" '$3 >= 35 && $15 >= 70' ${out} > ${out}_parsed

module purge && module load pigz
pigz -p8 ${out}
pigz -p8 ${out}_parsed
```

JobID: 50983004
