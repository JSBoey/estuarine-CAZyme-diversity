# 0. Set up

## 0a. Directories

**Main working directory**: `/nesi/nobackup/uoa00348/boey/2023-estuarine-cazyme-diversity/`

| Item | Location | Comment |
| :--- | :--- | :--- |
| MAGs | `/nesi/project/ga02676/Waiwera_project/3.Binning/2.Final_bins` | Original set of final bins curated by Dave, Sze and Carmen. Also contains CheckM results. |
| KoFamScan<br>(version 1.3.0) | `/nesi/project/uoa02469/Software/kofam_scan_v1.3.0` | Annie updated the database on 29 Mar 2023. |
| Transporter Classification Database<br>(TCDB; 12 Apr 2023) | `/nesi/project/uoa02469/Databases/TCDB_20230412` | Boey updated database on 13 Apr 2023. |
| dbCAN<br>(version 11) | `/nesi/project/uoa02469/Databases/dbCAN2_v11` | |
| dbCAN-sub<br>(downloaded 11 Aug 2022) | `/nesi/project/uoa02469/Databases/dbCAN-sub_20220811` | |
| CAZy<br>(10 Aug 2022) | `/nesi/project/uoa02469/Databases/CAZyDB_20220806` | Compiled by Yin et al. of dbCAN. |
| SulfAtlas<br>(version 2.3.1) | `/nesi/project/uoa02469/Databases/SulfAtlas_v2.3.1` | Amino acid sequence database; HMM searchable [here](https://sulfatlas.sb-roscoff.fr/sulfatlashmm/). |
| InterProScan (version 5.61-93.0) | `bin/interproscan-5.61-93.0` | Module Java/17 needs to be loaded prior. Also, not all analyses ran, but seems to work fine for Pfam, TIGRFAM, and CDD. |

A soft link to main project directory was also set up.

```bash
ln -s /nesi/project/ga02676/Waiwera_project/boey_work/estuarine-CAZyme-diversity/ main
```

## 0b. Additional software

### InterProScan update

NeSI has InterProScan as a module, but their version is 3 years old. Best to have a newer version. 

These instructions were adapted from [here](https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html).


```bash
wget --directory-prefix bin/ https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.61-93.0/interproscan-5.61-93.0-64-bit.tar.gz
wget --directory-prefix bin/ https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.61-93.0/interproscan-5.61-93.0-64-bit.tar.gz.md5

# Recommended checksum to confirm the download was successful:
md5sum -c interproscan-5.61-93.0-64-bit.tar.gz.md5
# Must return *interproscan-5.61-93.0-64-bit.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy.
```

\[2023-04-13 10:43\] `wget` is being run on a `screen` named `dl-ips`.

\[2023-04-14 14:26\] InterProScan version 5.61-93.0 installed in `bin`

### SignalP-6 GPU

Need to ask Dini to help with this. Otherwise, can continue to use CPU version, but only stick to CAZymes. Therefore, filter dbCAN identified CAZymes then run it through SignalP.

### DeepTMHMM

Potentially need to run this on their [cloud server](https://dtu.biolib.com/DeepTMHMM). This is not ideal, but I don't think I have a choice.

## 0c. Additional comments

Make sure to copy the notebook and `scripts/` to the backed up project directory at the end of the day.

I'll add a script to automate that.

Content of `scripts/0-backup_project.sh`:

```bash
#!/bin/bash -e

TARGET=/nesi/project/ga02676/Waiwera_project/boey_work/2023-estuarine-cazyme-diversity/

cp -v 2023-estuarine-cazyme-diversity.ipynb $TARGET
cp -v -r scripts/ $TARGET

```

Once there, make sure to commit to github from the project directory.

```bash
# Stage
git add --all

# Commit
git commit -m "Today's update"

# Push
git push
```

If a mistake was made:

```bash
# Reset (change HEAD back to good commit)
git reset --hard good_commit

# Push
git push -f origin <last_good_commit>:<branch>
```

## On `data/` and `results/`

`data/` is only selectively and manually backed-up to prevent bloating. Test resulting compressed filesize each time a data file needs to be backed-up. `data/` is never pushed to the git repository. `results/` is backed-up but not pushed to the git repository.

## 0. Sequence statistics

It may be important for down the line. Obtain length statistics for all analysed sequences.

```bash
module purge
module load SeqKit

# Protein predictions
seqkit fx2tab -n -l data/2.orf_prediction/allbins_pred.faa \
  | gzip -c \
  > results/allbins_pred.faa.length.gz

# Database sequences (from FASTA)
for i in data/0.db/*.faa; do
  base=$(basename $i .faa | sed -E 's/([A-Za-z]+).*/\1/g' )
  seqkit fx2tab -n -l $i | gzip -c > results/${base}.faa.length.gz
done
```

# 1. Get high quality bins

In the MAGs directory, use CheckM stats to filter to copy MAGs with $\gt$ 70% completeness and $\lt$ 5% contamination. These bins are copied into `data/0.bins/`.

```bash
bash scripts/1.get_bins.sh
```

Contents of `1-get_bins.sh`:

```bash
#!/bin/bash -e

# Subset MAGs based on CheckM stats

# Directories
DIRIN=/nesi/project/ga02676/Waiwera_project/3.Binning/2.Final_bins
DIROUT=data/1.bins

mkdir -p $DIROUT

# Variables
COMP=70 # Completeness
CONT=5  # Contamination

# Subset CheckM results and copy MAGs to data/
awk -v comp="$COMP" -v cont="$CONT" -F "\t" '($3 > comp) && ($4 < cont) {print $6}' $DIRIN/checkm_data.txt \
  | xargs -I {} cp $DIRIN/{}.fna $DIROUT/
```

Also get the CheckM outputs for these high quality bins.

```bash
DIRIN=/nesi/project/ga02676/Waiwera_project/3.Binning/2.Final_bins
DIROUT=data/1.bins

head -n 1 $DIRIN/checkm_data.txt > results/hq_bins.checkm_data.txt
for i in $DIROUT/Ww*.fna; do
  bin=$(basename ${i} .fna)
  grep "${bin}" $DIRIN/checkm_data.txt >> results/hq_bins.checkm_data.txt
done
```

# 2. Predict ORFs

Predict ORFs from MAGs using `prodigal` in an array.

```bash
sbatch scripts/2-ORF_prediction.sl
```

Contents of `2-ORF_prediction.sl`:

```bash
#!/bin/bash -e
#SBATCH --job-name=prodigal
#SBATCH --account=uoa00348
#SBATCH --time=2:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-250
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

# Reference
printf "%s was called on %s\n" "$(basename $0)" "$(date)"

# Modules
module purge
module load prodigal/2.6.3-GCCcore-7.4.0

# Directories
DIRIN=data/1.bins
DIROUT=data/2.orf_prediction

mkdir -p $DIROUT

# Variables
ARR=($DIRIN/*.fna)
INPUT=${ARR[$SLURM_ARRAY_TASK_ID]}
NAME=$(basename $INPUT .fna)
OUTPUT=${DIROUT}/${NAME}_pred

# Run prodigal
prodigal \
  -i ${INPUT} \
  -f gff \
  -a ${OUTPUT}.faa \
  -d ${OUTPUT}.fna \
  -o ${OUTPUT}.gff \
  -p single
```

**Job ID:** 34268687

ORF predictions were archived, compressed, and then backed-up.

```bash
module purge
module load pigz

mkdir -p main/data

cd data/2.orf_prediction

tar -cvf - Ww* | pigz -c -p 6 > ../../main/data/orf_prediction.tar.gz

cd ../../
```

## 2.1 Combine and clean up predictions

Concatenate all ORF predictions then remove metadata from the ORF predictions.

```bash
bash scripts/2.1-ORF_cleanup.sh
```

Contents of `2.1-ORF_cleanup.sh`:

```bash
#!/bin/bash -e

# Concatenate all ORF predictions and remove metadata from ORF predicted sequences.
# Metadata is also consolidated from the GFF files.

# Directories
DIR=data/2.orf_prediction

# Concatenate predictions and remove metadata
cat $DIR/*.faa \
  | cut -f 1 -d ' ' \
  > $DIR/allbins_pred.faa

# Create another without asterisk for interproscan
sed -e 's/\*//g' $DIR/allbins_pred.faa > $DIR/allbins_pred.noast.faa

# Consolidate GFF data
printf "bin\tnode\tsource\ttype\tstart\tend\tgff_score\tstrand\tphase\tseqid\tpartial\tstart_type\trbs_motif\trbs_spacer\tgc_cont\tconf\tscore\tcscore\tsscore\trscore\tuscore\ttscore\n" > $DIR/allbins_pred.metadata.tsv

for i in $DIR/*.gff; do
  bin=$(basename $i .gff)
  grep -v '#' $i \
    | sed -e 's/;$//' \
    | sed -e 's/;/\t/g' \
    | sed -E 's/\w+=//g' \
    | awk '{FS="\t"; OFS="\t"} {split($9, a, "_"); $1=$1"_"a[2]; print}' \
    | awk -v mag="$bin" '{FS="\t"; OFS="\t"} {print mag, $0}' \
    >> $DIR/allbins_pred.metadata.tsv
done
```

# 3. ORF annotation

Predict annotations for predicted ORFs against:
* KEGG (via KoFamScan)
* dbCAN
* CAZy
* SulfAtlas
* TCDB
* InterPro

Also predict signal peptides (via SignalP-6) and transmembrane proteins (via DeepTMHMM).

Keep in mind that with InterProScan, sequences will need to be chunked to around 80,000 sequences per file.

## 3.1 InterProScan

Need to split inputs into 80,000 sequences per file.

```bash
ml purge
ml SeqKit/2.2.0

mkdir -p data/tmp

seqkit split \
  data/2.orf_prediction/allbins_pred.faa \
  --out-dir data/tmp \
  --by-size 80000
```

Split into 10 files.

Proceed to InterPro annotations.

Contents of `3.1-interproscan.sl`:

```bash
#!/bin/bash -e
#SBATCH --job-name=ipr5
#SBATCH --account=uoa00348
#SBATCH --time=2:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-9
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

# Modules
module purge
module load Java/17

# Directories
INDIR=data/tmp
OUTDIR=data/3.annotation

mkdir -p ${OUTDIR}

# Variables
interproscan=bin/interproscan-5.61-93.0/interproscan.sh

ARR=(${INDIR}/allbins_pred.part*.faa)
INFILE=${ARR[$SLURM_ARRAY_TASK_ID]}
FILEPART=$(basename ${INFILE} .faa | sed -E 's/.*(part_[0-9]*)/\1/g')
OUTBASE=${OUTDIR}/allbins_pred.interpro.${FILEPART}
APPL=Pfam,TIGRFAM,CDD
FORMAT=xml,tsv
CPU=$(($SLURM_CPUS_PER_TASK - 2))


# Run InterProScan
$interproscan \
  --applications ${APPL} \
  --cpu ${CPU} \
  --formats ${FORMAT} \
  --input ${INFILE} \
  --output-file-base ${OUTBASE}
```

**Job ID:** 34314370

### Back up results

XML files are compressed then copied to project directory.

TSV files are concatenated, compressed, then copied to project directory.

```bash
ml purge
ml pigz

# Concatenate then compress TSV files
cat data/3.annotation/*.tsv | pigz -c -p 4 > results/allbins_pred.interpro.tsv.gz

# Archive then compress XML files
tar -cvf - data/3.annotation/*interpro.*.xml | pigz -p 4 > results/allbins_pred.interpro.xml.tar.gz
```

## 3.2 dbCAN

This will annotate against dbCAN and dbCAN-sub using HMMER3.

Contents of `3.2-dbcan.sl`:

```bash
#!/bin/bash -e
#SBATCH --job-name=dbcan
#SBATCH --account=uoa00348
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

# Modules
module purge
module load HMMER/3.3.2-GCC-11.3.0

# Directories
INDIR=data/2.orf_prediction
OUTDIR=data/3.annotation
DBDIR=/nesi/project/uoa02469/Databases

# Variables 
INFILE=${INDIR}/allbins_pred.faa
INBASE=$(basename ${INFILE} .faa)

## Database for dbCAN2
DB1=${DBDIR}/dbCAN2_v11/dbCAN-HMMdb-V11
## Output for dbCAN2
OUTFILE1=${OUTDIR}/${INBASE}.dbcan.domtbl

## Database for dbCAN-sub
DB2=${DBDIR}/dbCAN-sub_20220811/dbCAN_sub.hmm
## Output for dbCAN-sub
OUTFILE2=${OUTDIR}/${INBASE}.dbcan-sub.domtbl

## dbCAN2 size 
SZDB1=$(hmmstat ${DB1} | tail -n 1 | cut -f 1 -d ' ')

## dbCAN-sub size
SZDB2=$(hmmstat ${DB2} | tail -n 1 | cut -f 1 -d ' ')

# Run
## dbCAN
printf "[%s]\tStart dbCAN annotation. Database size = %s\n" "$(date)" $SZDB1

hmmsearch \
  -Z $SZDB1 \
  --cpu ${SLURM_CPUS_PER_TASK} \
  --domtblout ${OUTFILE1} \
  -o /dev/null \
  $DB1 \
  $INFILE

printf "[%s]\tFinished dbCAN annotation. Domain table saved to %s\n" "$(date)" ${OUTFILE1}

## dbCAN-sub
printf "[%s]\tStart dbCAN-sub annotation. Database size = %s\n" "$(date)" $SZDB2

hmmsearch \
  -Z $SZDB2 \
  --cpu ${SLURM_CPUS_PER_TASK} \
  --domtblout ${OUTFILE2} \
  -o /dev/null \
  $DB2 \
  $INFILE

printf "[%s]\tFinished dbCAN-sub annotation. Domain table saved to %s\n" "$(date)" ${OUTFILE2}

```

**Job ID:** 34412811

### Filter and backup

All dbCAN outputs were filtered using `hmmsearch-parser.sh`.

```bash
module purge
module load pigz

# Filter outputs
for i in data/3.annotation/*dbcan*; do
  hits=$(grep -c -v '#' ${i} | sed -E 's/.*:([0-9]+)/\1/g')
  printf "Filtering %s containing %s hits\n" "${i}" "${hits}"
  base=$(basename ${i} .domtbl)
  bin/hmmsearch-parser.sh $i > results/${base}_parsed.tsv
  fhits=$(wc -l results/${base}_parsed.tsv | cut -f 1 -d ' ')
  printf "Complete. Retained %s hits\n" "${fhits}"
  printf "Compressing results\n"
  pigz --best -v -p 6 results/${base}_parsed.tsv
done

# Backup outputs
for i in data/3.annotation/*dbcan*; do
  base=$(basename ${i})
  pigz --best -v -c -p 6 ${i} > main/data/${base}.gz
done
```

## 3.3 KEGG (KoFamScan)

Remember to funnel the temporary files to `data/tmp/`, then have `KoFamScan` delete the temporary files.

Contents of `3.3-kofamscan.sl`:

```bash
#!/bin/bash -e
#SBATCH --job-name=kofam
#SBATCH --account=uoa00348
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

# Modules
module purge
module load \
  Ruby/3.0.1-GCC-11.3.0 \
  HMMER/3.3.2-GCC-11.3.0 \
  Parallel/20220922

# Directories
INDIR=data/2.orf_prediction
OUTDIR=data/3.annotation
SOFTDIR=/nesi/project/uoa02469/Software/kofam_scan_v1.3.0
TMPDIR=data/tmp

# Variables
INFILE=${INDIR}/allbins_pred.faa
INBASE=$(basename ${INFILE} .faa)
OUTFILE=${OUTDIR}/${INBASE}.kofam.tsv

## Software variables
kofamscan=${SOFTDIR}/bin/exec_annotation

KOLIST=${SOFTDIR}/db/ko_list
PROFILE=${SOFTDIR}/db/profiles
FORMAT=detail-tsv

# Run
$kofamscan \
  --format=$FORMAT \
  --profile=$PROFILE \
  --ko-list=$KOLIST \
  --cpu=${SLURM_CPUS_PER_TASK} \
  --tmp-dir=$TMPDIR \
  ${INFILE} \
  > ${OUTFILE}

```

**Job ID:** 34372470

### Filter and backup

KoFamScan outputs were filtered to retain only significant hits as results. These results were then compressed and backed-up.

```bash
module purge
module load pigz

grep "\*" data/3.annotation/allbins_pred.kofam.tsv \
  | pigz -c -p 4 \
  > results/allbins_pred.kofam_sig.tsv.gz
```

KoFamScan outputs were also compressed and backed-up in the main directory under `data/`.

```bash
module load pigz

pigz --best -c -p 6 data/3.annotation/allbins_pred.kofam.tsv > main/data/allbins_pred.kofam.tsv.gz
```

## 3.4 Sequence alignment: Transporters, sulfatases, peptidases, and CAZy

Use DIAMOND `blastp` to detect:

- Transporters (TDCB)
- Sulfatases (SulfAtlas)
- Peptidases (MEROPS) \[The FTP site at EBI's end is down, cannot download database @ 17 Apr 2023\]
- CAZymes (CAZyDB)

### Copy and index databases in scratch space

```bash
# Copy existing database
mkdir data/0.db

cp /nesi/project/uoa02469/Databases/{CAZy,SulfAtlas,TCDB_2023}*/*.faa data/0.db

# Use only nometa for SulfAtlas
rm data/0.db/sulfatlas_v2.3.1.faa

# Index databases
module purge
module load DIAMOND/2.1.1-GCC-11.3.0

cd data/0.db

for i in *.faa; do
  dbbase=$(basename ${i} .faa)
  diamond makedb --in ${i} --db ${dbbase}.dmnd --threads 6 --verbose
done
```

Contents of `3.4-diamond.sl`:

```bash
#!/bin/bash -e
#SBATCH --job-name=dmnd
#SBATCH --account=uoa00348
#SBATCH --time=3:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-2
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

# Modules
module purge
module load DIAMOND/2.1.1-GCC-11.3.0

# Directories
INDIR=data/2.orf_prediction
OUTDIR=data/3.annotation
DBDIR=data/0.db

# Variables
DBARRAY=(${DBDIR}/*.dmnd)
DB=${DBARRAY[$SLURM_ARRAY_TASK_ID]}
DBBASE=$(basename ${DB} | sed -E 's/([A-Za-z]+).*/\1/g')

INFILE=${INDIR}/allbins_pred.faa
INBASE=$(basename "${INFILE}" .faa)
OUTFILE=${OUTDIR}/${INBASE}.${DBBASE}.tsv

## Software variables
FORMAT='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp ppos'
MAXTARGET=0

# Run
printf "[%s] diamond blastp against %s\n" "$(date)" "${DBBASE}"

diamond blastp \
  --very-sensitive \
  --threads ${SLURM_CPUS_PER_TASK} \
  --db ${DB} \
  --query ${INFILE} \
  --max-target-seqs ${MAXTARGET} \
  --outfmt ${FORMAT} \
  --out ${OUTFILE}

```

**Job ID:** 34434621

### Filter and backup

All sequence alignment outputs were filtered to retain hits that fulfilled the following criteria:

- Percent identity $\ge$ 30%
- Query coverage $\ge$ 70%
- E-value $\lt$ 0.001

All software outputs and filtered hits were backed up to `main/`.

```bash
module purge
module load pigz

# Filter outputs
for i in data/3.annotation/*{CAZy,sulf,tcdb}*; do
  hits=$(wc -l ${i} | cut -f 1 -d ' ')
  printf "Filtering %s containing %s hits\n" "${i}" "${hits}"
  base=$(basename ${i} .tsv)
  awk -F '\t' -v id=30 -v qcov=70 -v evalue=0.001 \
    '$3 >= id && $11 < evalue && $13 > qcov' \
    ${i} \
    > results/${base}_filt.tsv
  fhits=$(wc -l results/${base}_filt.tsv)
  printf "Complete. Retained %s hits\n" "${fhits}"
  printf "Compressing results\n"
  pigz --best -v -p 6 results/${base}_filt.tsv
done

# Backup outputs
for i in data/3.annotation/*{CAZy,sulf,tcdb}*; do
  base=$(basename ${i})
  pigz --best -v -c -p 6 ${i} > main/data/${base}.gz
done
```

## 3.5 Signal peptides

This is achieved using SignalP6.

```bash
#!/bin/bash -e
#SBATCH --job-name=signalp6-gpu
#SBATCH --account=uoa00348
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=12
#SBATCH --gpus-per-node=1
#SBATCH --output=slurm_out/%x.%j.out
#SBATCH --error=slurm_err/%x.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz

# Modules
module purge
module load \
  SeqKit/2.2.0 \
  Python/3.9.9-gimkl-2020a \
  CUDA/12.0.0

# Environment
source /nesi/project/uoa02469/Software/signalp6_fast-gpu/bin/activate

# Directories
INDIR=data/2.orf_prediction
OUTDIR=data/3.annotation/signalp6
TMPDIR=tmp/signalp6_$SLURM_JOB_ID

mkdir -p {$OUTDIR,$TMPDIR}

# Variables
PREDFILE=$INDIR/allbins_pred.faa
FORMAT=none
ORGANISM=other
MODE=fast
BATCH_SIZE=500

# Run
## Split input into multiple files
seqkit split2 -s 50000 $PREDFILE -O $TMPDIR

## Run SignalP6 on part files
for i in $TMPDIR/*.faa; do

  # Get part name
  part=$(basename $i .faa)

  # Create part result directory
  mkdir -p $TMPDIR/$part

  # Run
  signalp6 \
    --fastafile $i \
    --output_dir $TMPDIR/$part \
    --format $FORMAT \
    --organism $ORGANISM \
    --mode $MODE \
    --bsize $BATCH_SIZE \
    --write_procs $SLURM_CPUS_PER_TASK

  # Add part as prefix
  for j in $TMPDIR/$part/*; do

    # Get result basename
    filename=$(basename $j)

    # Move results to OUTDIR with part prefix
    mv $j $OUTDIR/${part}_${filename}

  done

done

```

**Job ID**: 35006571

Concatenate human readable outputs to `results/`

```bash
SIGNALP_DIR=data/3.annotation/signalp6

# Set-up output files
for i in $SIGNALP_DIR/allbins_pred.part_001*.{txt,gff3}; do
  suffix=$(basename $i | sed -E 's/allbins_pred.part_[0-9]+_//g')
  grep "#" $i > results/allbins_pred.signalp_${suffix}
done

# Concatenate 
for i in $SIGNALP_DIR/allbins_pred.*.{txt,gff3}; do
  suffix=$(basename $i | sed -E 's/allbins_pred.part_[0-9]+_//g')
  prefix=$(echo $i | sed -e "s/${suffix}//g")
  # Remove headers and "OTHER" from prediction_results.txt
  grep -v "#" ${prefix}${suffix} \
    | grep -v 'OTHER' \
    >> results/allbins_pred.signalp_${suffix}
done

```

### Backup

```bash
cd results/

module purge
module load pigz
pigz --best -p 8 *.{txt,gff3}

cd ../data/3.annotation/signalp6
tar -cvf - * | pigz --best -p 8 -c > ../../../main/data/allbins_pred.signalp.tar.gz
```

# 4. Gene and transcript count

Using Bowtie2 to map metagenomic and metatranscriptomic reads back to bins. Then use Subread's featureCounts to generate counts.

Create directory for this section.

```bash
mkdir -p data/4.feature_count/{index,bam,count}
mkdir -p data/0.w{g,t}s
```

## 4.1. Index scaffolds

```bash
# Load module
module purge
module load Bowtie2/2.4.5-GCC-11.3.0

# Create comma-separated line of fasta files
files=$(find data/1.bins -type f -printf '%p,' | sed 's/,$//')

# Index bins
bowtie2-build --threads 8 $files data/4.feature_count/index/allbins_scaffolds
```

## 4.2. Map gene and transcript reads

Copy the relevant read files from project directory:

- Whole genome shotgun: `/nesi/project/ga02676/Waiwera_project/1.Qual_filtered_trimmomatic`
- Whole transcriptome shotgun: `/nesi/project/ga02676/Waiwera_project/sze_works/1.WTS/3.trimmomatic`

```bash
rsync -av /nesi/project/ga02676/Waiwera_project/1.Qual_filtered_trimmomatic/*_R{1,2}* data/0.wgs
rsync -av /nesi/project/ga02676/Waiwera_project/sze_works/1.WTS/3.trimmomatic/*.fastq.gz data/0.wts
```

The file names need to be standardised due to WTS and WGS having different naming conventions. This is the convention to follow:

```
<read_type>.<sample_type>.<site>_<biological_replicate>.<read_pair>.fastq.gz
```

Furthermore, water samples in WGS reads were run on separate lanes. Remember to concatenate files from the same samples.

```bash
# Working on WGS reads
cd data/0.wgs

# Concatenate reads from same water samples
module purge
module load pigz

for i in {1..9}; do
  for j in {1..2}; do
    cat Filt.S${i}_L?_R${j}.fastq.gz > Filt.S${i}_R${j}.fastq.gz
  done
done

# Remove individual lane reads
rm Filt.S?_L?_R?.fastq.gz

# Rename
for i in Filt*; do
  read_type=WGS
  sample_type=Filt
  site=$(echo $i | sed -E 's/.*S([0-9]).*/\1/g')
  biol_rep=1
  read_pair=$(echo $i | sed -E 's/.*(R[0-9]).*/\1/g')
  newname=${read_type}.${sample_type}.S${site}_${biol_rep}.${read_pair}.fastq.gz
  mv $i $newname
done

for i in Sed*; do
  read_type=WGS
  sample_type=Sed
  site=$(echo $i | sed -E 's/.*S([0-9]).*/\1/g')
  biol_rep=$(echo $i | sed -E 's/.*Sample([0-9]).*/\1/g')
  read_pair=$(echo $i | sed -E 's/.*(R[0-9]).*/\1/g')
  newname=${read_type}.${sample_type}.S${site}_${biol_rep}.${read_pair}.fastq.gz
  mv $i $newname
done

# Working on WTS reads
cd ../0.wts

# Rename
for i in *Filt*; do
  read_type=WTS
  sample_type=Filt
  site=$(echo $i | sed -E 's/.*Filt([0-9]).*/\1/g')
  biol_rep=1
  read_pair=$(echo $i | sed -E 's/.*(R[0-9]).*/\1/g')
  newname=${read_type}.${sample_type}.S${site}_${biol_rep}.${read_pair}.fastq.gz
  mv $i $newname
done

for i in *Sed*; do
  read_type=WTS
  sample_type=Sed
  site=$(echo $i | sed -E 's/.*S([0-9]).*/\1/g')
  biol_rep=$(echo $i | sed -E 's/.*Sed([0-9]).*/\1/g')
  read_pair=$(echo $i | sed -E 's/.*(R[0-9]).*/\1/g')
  newname=${read_type}.${sample_type}.S${site}_${biol_rep}.${read_pair}.fastq.gz
  mv $i $newname
done

```

Contents of `scripts/4.2-map_reads.sl`:

```bash
#!/bin/bash -e
#SBATCH --job-name=bt2map
#SBATCH --account=uoa00348
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=36
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-65
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz
#SBATCH --partition=milan

# Modules
module purge
module load \
  Bowtie2/2.4.5-GCC-11.3.0 \
  SAMtools/1.16.1-GCC-11.3.0

# Directories
WGSDIR=data/0.wgs
WTSDIR=data/0.wts
INDEXDIR=data/4.feature_count/index
BAMDIR=data/4.feature_count/bam

export TMPDIR=data/tmp/bt2map/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $TMPDIR

# Variables
ARR=(${WGSDIR}/*.R1.fastq.gz ${WTSDIR}/*.R1.fastq.gz)
INBASE=$(basename ${ARR[$SLURM_ARRAY_TASK_ID]} .R1.fastq.gz)
INDIR=$(dirname ${ARR[$SLURM_ARRAY_TASK_ID]})
INDEX=${INDEXDIR}/allbins_scaffolds
READ1=${INDIR}/${INBASE}.R1.fastq.gz
READ2=${INDIR}/${INBASE}.R2.fastq.gz
OUTBAM=${BAMDIR}/${INBASE}.bam

# Run
bowtie2 -p $SLURM_CPUS_PER_TASK -x $INDEX -1 $READ1 -2 $READ2 \
  | samtools view -@ $SLURM_CPUS_PER_TASK -bS - \
  | samtools sort -@ $SLURM_CPUS_PER_TASK -o $OUTBAM -

# Clean up
rm -rf $TMPDIR

```

**Job ID**: 34936855

## 4.3. Count features

Use `Subread`'s `featureCounts` for this. Start by prepare simple annotation format (SAF) from compiled GFF.

```bash
# Load modules
module purge
module load Subread/2.0.3-GCC-11.3.0
module load pigz

# Directories
ORFDIR=data/2.orf_prediction
COUNTDIR=data/4.feature_count/count
BAMDIR=data/4.feature_count/bam

# Create SAF
printf "GeneID\tChr\tStart\tEnd\tStrand\n" > $COUNTDIR/allbins_pred.saf

cat $ORFDIR/allbins_pred.metadata.tsv \
  | tail -n +2 \
  | cut -f 2,5,6,8 \
  | awk '{FS="\t"; OFS="\t"} {print $1, $gsub("_[0-9]+$", "", $1), $2, $3, $4}' \
  >> $COUNTDIR/allbins_pred.saf

# Count genes and transcripts
for i in WGS WTS; do
  featureCounts -p -T 12 \
    -a $COUNTDIR/allbins_pred.saf \
    -F SAF \
    -o $COUNTDIR/${i}_count.tsv \
    $BAMDIR/${i}.*.bam
done

# Clean up counts and funnel to results/
for i in WGS WTS; do
  grep -v '#' $COUNTDIR/${i}_count.tsv \
    | sed -e 's/data\/4\.feature_count\/bam\///g' \
    | sed -e 's/\.bam//g' \
    | sed -e "s/${i}\.//g" \
    | pigz -p 8 --best -c \
    > results/${i}_clean_count.tsv.gz
done
```

## Backup

Need to backup all aligned reads and count files. BAM files are still quite big, so I will try to create CRAM files for backup purposes.

### Back-up aligned reads

```bash
# Set-up directories and variables
mkdir data/4.feature_count/cram

BAMDIR=data/4.feature_count/bam
CRAMDIR=data/4.feature_count/cram

# Load SAMtools
module purge
module load SAMtools/1.16.1-GCC-11.3.0

# Create the reference
cat data/1.bins/Ww*.fna > $CRAMDIR/allbins_scaffolds.fna
cd $CRAMDIR

samtools faidx allbins_scaffolds.fna

cd ../../../
```

Conversion of BAM to CRAM is performed using `scripts/4-bam2cram.sl`:

```bash
#!/bin/bash -e
#SBATCH --job-name=bam2cram
#SBATCH --account=uoa00348
#SBATCH --time=30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-65
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jian.sheng.boey@auckland.ac.nz
#SBATCH --partition=milan

# Modules
module purge
module load SAMtools/1.16.1-GCC-11.3.0

# Directories
BAMDIR=data/4.feature_count/bam
CRAMDIR=data/4.feature_count/cram

# Variables
ARR=($BAMDIR/*.bam)
INFILE=${ARR[$SLURM_ARRAY_TASK_ID]}
INBASE=$(basename ${INFILE} .bam)
OUTFILE=${CRAMDIR}/${INBASE}.cram

# Run
samtools view -C \
  -@ $SLURM_CPUS_PER_TASK \
  -T ${CRAMDIR}/allbins_scaffolds.fna \
  -o ${OUTFILE} \
  ${INFILE}

```

**Job ID**: 34958905

> Remember to convert CRAM back to BAM in order to regenerate count files if necessary.
> The reference seqeuences (compressed) are stored with the CRAM files.

Archive files to `main/`

```bash
cd data/4.feature_count/cram

# Compress reference
module load pigz
pigz --best -p 8 allbins_scaffolds.fna

tar -cvf - * > ../../../main/data/read_alignment_cram.tar
```

### Back-up count files and summaries

```bash
cd data/4.feature_count/count

module load pigz
tar -cvf - * | pigz --best -c -p 8 > ../../../main/data/gene_transcript_counts.tar.gz
```

# 5. Read coverage

This was done after the data was generated. Therefore, the data was pulled from the project directory.

## 5.1 Create directories and pull inputs

Create a separate directory in a scratch directory of choice.

```sh
# Create necessary directories
mkdir -p {alignments/{bam,cram},coverage,counts}

# Pull inputs
main=/nesi/project/ga02676/Waiwera_project/boey_work/estuarine-CAZyme-diversity
cp ${main}/data/{read_alignment_cram.tar,gene_transcript_counts.tar.gz} .

# Decompress and extract
tar -xvf read_alignment_cram.tar -C alignments/cram
tar -xzvf gene_transcript_counts.tar.gz -C coverage/
```

`gene_transcript_counts.tar.gz` is required as it contains the original SAF file used in featureCounts. This is the basis of the BED file.

## 5.2 File conversions

BEDTools does not accept CRAM inputs. Therefore all CRAM files must be converted to BAM files as per the following file.

`cram2bam.sl` as follows:

```sh
#!/bin/bash -e
#SBATCH --job-name=cram2bam
#SBATCH --account=ga02676
#SBATCH --time=30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-65

module purge
module load SAMtools/1.19-GCC-12.3.0

arr=(alignments/cram/*.cram)
in=${arr[$SLURM_ARRAY_TASK_ID]}
bn=$(basename ${in} .cram)
outdir=alignments/bam
ref=$(dirname ${in})/allbins_scaffolds

printf "Input file: %s\n" ${in}

samtools view \
  -t ${ref}.fna.fai -T ${ref}.fna \
  -@ $SLURM_CPUS_PER_TASK \
  -b -o ${outdir}/${bn}.bam \
  ${in}

printf "Output file: %s\n" ${outdir}/${bn}.bam
```

SAF file must be converted to BED file. However, the ordering needs to be in the order of the BAM files as genes/meta-features were ordered based on the indices of the metagenomic assemblies. Thus, coordinate and name ordering is not relevant here. **This step is important to use `bedtools coverage -sorted` for memory efficiencies and speed-up.**

`saf2bed.sh` as follows:

```sh
#!/bin/bash -e

module purge
module load SAMtools/1.19-GCC-12.3.0 BEDTools/2.30.0-GCC-11.3.0

bn=$(basename $1 .saf)
dn=$(dirname $1)

# Convert SAF to BED
sed '1d' $1 \
  | awk 'BEGIN {FS="\t"; OFS="\t"} {print $2, $3, $4, $1, $5}' \
  > ${dn}/${bn}.bed

# Compare 2 random BAM to ensure similar NODE ordering
sam1=$(ls alignments/bam/WGS.*.bam | shuf -n 1)
sam2=$(ls alignments/bam/WTS.*.bam | shuf -n 1)

samtools view -H ${sam1} | grep '@SQ' | cut -f 2 | sed 's/SN://g' > names1
samtools view -H ${sam2} | grep '@SQ' | cut -f 2 | sed 's/SN://g' > names2

# Sort BED based on BAM order
if cmp -s names1 names2; then
  # BEDTools require names.txt with columns: chr length
  samtools view -H ${sam1} \
    | grep '@SQ' \
    | cut -f2,3 \
    | sed 's/[SL]N://g' \
    > names.txt
  bedtools sort -g names.txt -i ${dn}/${bn}.bed > ${dn}/${bn}.sorted.bed
  rm names{1,2}
else
  printf "Nodes are sorted differently"
fi

# Sort BED based on chromosome names and start positions (not run)
# sort -k1,1 -k2,2n ${dn}/${bn}.bed > ${dn}/${bn}.sorted.bed

module purge
unset sam1 sam2
```

## 5.3 Calculate coverage per sample

`feature_cov.sl` as follows:

```sh
#!/bin/bash -e
#SBATCH --job-name=fcov
#SBATCH --account=ga02676
#SBATCH --time=30:00
#SBATCH --partition=milan
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-65

module purge
module load BEDTools/2.30.0-GCC-11.3.0

arr=(alignments/bam/*.bam)
in=${arr[$SLURM_ARRAY_TASK_ID]}
bn=$(basename ${in} .bam)

printf "Input: %s\n" ${in}

printf "Calculating coverage\n"

bedtools coverage \
  -sorted -g names.txt \
  -a counts/allbins_pred.sorted.bed \
  -b ${in} \
  > coverage/${bn}.cov.tsv
```

The script above obtains coverage per gene, for relative activity, I need coverage per contig as well:

`contig_cov.sl`

```sh
#!/bin/bash -e
#SBATCH --job-name=ccov
#SBATCH --account=ga02676
#SBATCH --time=30:00
#SBATCH --partition=milan
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-65

module purge
module load SAMtools/1.19-GCC-12.3.0

arr=(alignments/bam/*.bam)
in=${arr[$SLURM_ARRAY_TASK_ID]}
bn=$(basename ${in} .bam)

printf "Input: %s\n" ${in}

printf "Calculating coverage\n"

samtools coverage -o coverage/${bn}.contig_cov.tsv ${in}
```

Then get only relevant columns:

`clean_contig_cov.sh`

```sh
#!/bin/bash -e

xjoin() {
    local f

    if [ "$#" -lt 2 ]; then
        echo "xjoin: need at least 2 files" >&2
        return 1
    elif [ "$#" -lt 3 ]; then
        join "$1" "$2"
    else
        f=$1
        shift
        join "$f" <(xjoin "$@")
    fi
}

for i in coverage/*.contig_cov.tsv; do
  bn=$(basename $i .contig_cov.tsv)
  printf "Processing %s\n" ${bn}
  (head -n 1 ${i} && tail -n +2 ${i} | sort) \
    | cut -f 1,4,5 \
    | sed -e "s/numreads/${bn}.numreads/" -e "s/covbases/${bn}.covbases/" \
    | sed -e 's/ /\t/g' \
    > coverage/${bn}.contig_scov.tsv
done

for i in WGS WTS; do
  printf "Joining %s\n" ${i}
  xjoin coverage/${i}.*.contig_scov.tsv \
    | sed 's/#//g' \
    > ${i}.contig_coverage.tsv
done

printf "Completed\n"
```

## 5.4 Join columns

Outputs from `bedtools coverage` are the following columns:

1. Contig (Feature ala BEDTools)
2. ORF start
3. ORF end
4. Gene (Meta-feature ala BEDTools)
5. Strand
6. Number of aligned reads overlapping gene
7. Number of bases of all aligned reads covering gene (base coverage per gene)
8. Length of gene
9. Fraction of bases of all aligned reads covering gene (read coverage per gene)

Only columns 4, 6, and 7 are extracted for further use. Column headers are appended as per `clean_fcov.sh`:

```sh
#!/bin/bash -e

# Select relevant columns in output and append headers

arr=(coverage/*.cov.tsv)

for file in ${arr[@]}; do
  nm=$(basename ${file} .cov.tsv)
  printf "Appending %s\n" ${nm}
  cut -f 4,6,7 ${file} \
    | sed "1i\node\t${nm}.reads_overlapped\t${nm}.bases_covered" \
    > coverage/${nm}.scov.tsv
done
```

Then, a recursive join is used to obtain a table of read depth and bases covered via `recursive_join.sh`:

```bash
#!/bin/bash -e

# Define recursive function to join multiple scov tables
# Thanks to rudimeier from StackOverflow:
# https://unix.stackexchange.com/questions/364735/merge-multiple-files-with-join
xjoin() {
    local f

    if [ "$#" -lt 2 ]; then
        echo "xjoin: need at least 2 files" >&2
        return 1
    elif [ "$#" -lt 3 ]; then
        join "$1" "$2"
    else
        f=$1
        shift
        join "$f" <(xjoin "$@")
    fi
}

for i in WGS WTS; do
  xjoin coverage/${i}*.scov.tsv \
    | sed "s/${i}.//g" \
    | sed 's/ /\t/g' \
    > ${i}.feature_coverage.tsv

done
```

## 5.5 Backup and transfer

Relevant tables (`*.feature_coverage.tsv`) transferred to working local directory using `rclone`:

```sh
module load rclone

for i in *.feature_coverage.tsv; do
  rclone copy --verbose ${i} <destination>:<directory>
done
```

Original outputs archived, compressed, and transferred to project directory:

```bash
module load pigz
tar -vc coverage/*.cov.tsv coverage/*.contig_cov.tsv *.feature_coverage.tsv *.contig_coverage.tsv \
  | pigz -p 8 \
  > $main/data/gene_transcript_coverage.tar.gz
```

Scripts transferred to project directory

```bash
cp *.s{l,h} $main/scripts/
```

## 6. Read lengths

Required for read coverage calculations

Use BAM files generated prior

`calc_read_length.sl`

```sh
#!/bin/bash -e
#SBATCH --job-name=readlen
#SBATCH --account=ga02676
#SBATCH --time=20:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_out/%x.%j.%A_%a.out
#SBATCH --error=slurm_err/%x.%j.%A_%a.err
#SBATCH --array=0-65

module purge
module load SAMtools/1.19-GCC-12.3.0

mkdir -p alignments/read_length

arr=(alignments/bam/*.bam)
file=${arr[$SLURM_ARRAY_TASK_ID]}
base=$(basename ${file} .bam)
out=alignments/read_length/${base}.rlen.tsv

printf "Working on %s\n" ${file}
printf "read_length\tcount\n" > ${out}

samtools stats -@ $SLURM_CPUS_PER_TASK ${file} \
  | grep ^RL \
  | cut -f 2- \
  >> ${out}
```

Concatenate individual files into one big file with relevant sample information

`cat_read_length.sh`

```sh
#!/bin/bash -e

printf "srctype\tsample\tread_length\tcount\n" > read_lengths.csv

for file in alignments/read_length/*.tsv; do
  name=$(basename ${file} .rlen.tsv)
  srctype=$(echo ${name} | sed 's/\([A-Z]\+\)\..*/\1/g')
  sample=$(echo ${name} | sed 's/[A-Z]\+\.\(.*\)/\1/g')
  printf "Joining %s %s\n" ${srctype} ${sample}
  tail -n+2 ${file} \
    | sed "s/^/${srctype}\t${sample}\t/g" \
    >> read_lengths.csv
done
```

----


