# Simulate RNA-Seq Reads Using Polyester

Snakemake workflow to simulate RNA-Seq reads from fasta files using the Polyester R package.

To run, you'll need conda and snakemake > 5.10 installed. 

    If you don't have conda or snakemake:
        1. Download and install miniconda [here](https://docs.conda.io/en/latest/miniconda.html)
        2. Create a conda environment with snakemake
        """
        conda create -n snakemake-env snakemake==5.13
        """
        3. Activate that enviroment
        """
        conda activate snakemake-env
        """

    To test the workflow, run `snakemake -s simreads.snakefile --use-conda --configfiles config/test.yml`

    To run other datasets, build a configfile for your dataset of choice.


## Simulate from all transcripts in a fasta file (no gff file required)

    """
    snakemake -s simreads.snakefile --use-conda --configfiles config/QfO.yml -p
    """

## Simulate only from certain features, gff required
      
    """    
    snakemake -s simreads.snakefile --use-conda --configfiles config/Hsapiens_ensembl.yml -p
    """
    
    Note: this is currently limited to the `gene_name` attribute, but can be modified for additional attributes

