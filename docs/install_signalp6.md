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













Directories for CPU and GPU versions were made.

```bash
mkdir signalp6_fast-{c,g}pu
```

The directory containing all necessary install files of unpacked `signalp6_fast` were copied to `signalp6_fast-{c,g}pu`.

```bash
cp -r signalp6_fast/* signalp6_fast-cpu
cp -r signalp6_fast/* signalp6_fast-gpu
```
