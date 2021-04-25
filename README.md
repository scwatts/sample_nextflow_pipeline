# A simple nextflow pipeline
This repository contains a simple nextflow pipeline to demonstrate basic features of the framework and DSL2. This example is
intended to serve as a tutorial.


To run this pipeline, first install dependencies through conda and activate the environment:
```bash
conda create -p $(pwd -P)/conda_env/ --yes nextflow skesa biopython
conda activate conda_env/
```

Launch it 🚀:
```bash
./tutorial_pipeline.nf
```
