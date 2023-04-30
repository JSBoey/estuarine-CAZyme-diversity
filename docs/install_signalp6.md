# Install process for SignalP6 (CPU and GPU)

## 1. Set up
The package was downloaded into `bin/`, then copied to `/nesi/project/uoa02469/Software/` and unpacked.

```bash
cp bin/signalp-6.0g.fast.tar.gz /nesi/project/uoa02469/Software
cd /nesi/project/uoa02469/Software/
tar -xvzf signalp-6.0g.fast.tar.gz
```

Python 3.9.9 was loaded and then a virtual environment was created for CPU and GPU versions of the software.

```bash
module load Python/3.9.9-gimkl-2020a

python -m venv signalp6_fast-cpu
python -m venv signalp6_fast-gpu
```

## 2. Install CPU version

```bash
# Activate signalp6_fast-cpu environment
source signalp6_fast-cpu/bin/activate

# Install using pip
cd signalp6_fast
pip install signalp-6-package/

# Test by asking for help
signalp6 -h

# If all is good, the usage help should appear.

# Copy models into signalp
cp signalp-6-package/models/* ../signalp6_fast-cpu/lib/python3.9/site-packages/signalp/model_weights

# Check models are present
ls ../signalp6_fast-cpu/lib/python3.9/site-packages/signalp/model_weights
```

Usage for CPU version involves 

1. Activating the environment: `source /nesi/project/uoa02469/Software/signalp6_fast-cpu/bin/activate`
2. Calling signalp6 directly from terminal: `signalp6 -h`
3. When done, enter command: `deactivate`

### Testing

In `test/`, there are two fasta files:

- test_10.faa contains 9 sequences (don't ask me why)
- test_1k.faa contains 914 sequences (again, idk...)
- test.faa contains 36723 sequences

Running from the main working directory `/nesi/nobackup/uoa00348/boey/2023-estuarine-cazyme-diversity`:

```bash
signalp6 -ff test/test_10.faa -od test/ -f none
```

This produces intended outputs:

```
output.gff3  
output.json  
prediction_results.txt  
processed_entries.fasta  
region_output.gff3
```

## Install GPU version

Run the following to obtain a interactive shell session with a GPU node

```bash
srun --account uoa00348 --job-name "test-sp6gpu" --gpus-per-node A100-1g.5gb:1 --cpus-per-task 8 --mem 16GB --time 2:00:00 --pty bash
```

Load CUDA

```bash
module load CUDA/12.0.0
```

Install as per CPU version

```bash
cd /nesi/project/uoa02469/Software/

# Activate signalp6_fast-cpu environment
source signalp6_fast-gpu/bin/activate

# Install using pip
cd signalp6_fast
pip install signalp-6-package/

# Test by asking for help
signalp6 -h

# Copy models into signalp
cp signalp-6-package/models/* ../signalp6_fast-gpu/lib/python3.9/site-packages/signalp/model_weights

# Check models are present
ls ../signalp6_fast-gpu/lib/python3.9/site-packages/signalp/model_weights
```

Convert models

```bash
signalp6_convert_models gpu
```

## Compare SignalP CPU vs GPU

### Create test files

```bash
cd test/

module purge
module load SeqKit/2.2.0

for n in 100 5000 25000; do
  seqkit sample -n $n ../data/2.orf_prediction/allbins_pred.faa -o test_tmp.faa
  seqnum=$(grep -c '>' test_tmp.faa)
  mv test_tmp.faa test_seqsize_${seqnum}.faa
done
```

### Test GPU

Start a GPU interactive session.

```bash
srun \
  --account uoa00348 \
  --job-name "InteractiveGPU" \
  --gpus-per-node A100-1g.5gb:1 \
  --cpus-per-task 12 \
  --mem 16GB \
  --time 1:00:00 \
  --pty bash
```

Activate source for SignalP6-GPU

```bash
source /nesi/project/uoa02469/Software/signalp6_fast-gpu/bin/activate
```

Test per sequence file

```bash
module load CUDA/12.0.0 Python/3.9.9-gimkl-2020a

for file in test_seqsize_*.faa; do
  seqnum=$(echo $file | sed -E 's/.*_([0-9]+).*/\1/g')
  mkdir -p GPU_${seqnum}
  printf "Trying %s sequences\n" "${seqnum}"
  time signalp6 -ff $file -od GPU_${seqnum} -f none -bs 500 -wp 12
done
```

Results for GPU

| n<sub>seq</sub> | Time | Batch size |
| :-------------- | :--- | :--------- |
| 24841 | 12 minutes | 100 | 
| 5000  | 2.5 minutes | 100 |
| 87    | 8 seconds | 100 |   

GPU is doing about 20-40 sequences per second.

### Test CPU

Activate source for SignalP6-GPU

```bash
source /nesi/project/uoa02469/Software/signalp6_fast-cpu/bin/activate
```

Test per sequence file

```bash
for file in test_seqsize_{87,5000}.faa; do
  seqnum=$(echo $file | sed -E 's/.*_([0-9]+).*/\1/g')
  mkdir -p CPU_${seqnum}
  printf "Trying %s sequences\n" "${seqnum}"
  time signalp6 -ff $file -od CPU_${seqnum} -f none -bs 100 -wp 12
done
```

Results for CPU

| n<sub>seq</sub> | Time | Batch size |
| :-------------- | :--- | :--------- |
| 24841 | ~ 1.5 hours | 100 |
| 5000  | ~ 30 mins   | 100 |
| 87    | 23 seconds  | 100 |

In general, CPU is predicting at about 1-6 sequences per second. My guess is that it depends on the sequence itself.



