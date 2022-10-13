# metagenomic_anomaly_detector

A generalizable framework to create anomaly-detecting DeepSVDD models of metagenomic samples. takes .fasta/.fastq formatted inputs and returns a model and/or scores from the model on test samples. GPU acceleration will be automatically employed if available at training time.

The tool can be constructed using conda simply with the provided .yml.
* If you dont have miniconda, see installation here: https://docs.conda.io/en/latest/miniconda.html
    * Download the latest linux version
    * Navigate to the file you installed (should be named 'Miniconda3-latest-Linux-x86_64.sh',likely in the Downloads folder) and execute the following in terminal:
    ```
    ./Miniconda3-latest-Linux-x86_64.sh
    ```
Follow and agree to the installation prompts. You may need to restart your terminal session at the conclusion of installation.

Enter the base conda environment, and create the metagenomic anomaly detector conda environment from the .yml
```
conda activate base
conda env create -f metagenomic_anomaly_detector.yml
conda activate metagenomic_anomaly_detector
```

From here, construct the following expected folders for the command defaults (note: if you specify folders at runtime this is not necessary)
```
mkdir anomaly_dir
mkdir control_dir
mkdir models
mkdir test_samples
mkdir test_sample_scores
```

To complete set up, install the PanGIA database by running the following command:
```
wget PLACEHOLDER_LINK_TO_PANGIA_DB ./tools/pangia_db/
```

# The Workflow Overview

![](images/MBM_flowchart.png)

The main set of analysis is carried out by running three scripts in succession:
* run_pangia_jellyfish.py
* compile_and_compress_outputs.py
* train_DeepSVDD.py

Which will in full, from fasta/fastq inputs in specified folder locations, generate feature vectors and train DeepSVDD models over a small grid search to return an optimally fit model. A fourth optional script can be run with a trained model to evaluate testing data after this process is completed:
* get_score_from_model.py

Each script comes with command line arguments that can be checked by passing --help. The processes are detailed script-by-script below.

# run_pangia_jellyfish.py
# compile_and_compress_outputs.py
# train_DeepSVDD.py
# get_score_from_model.py