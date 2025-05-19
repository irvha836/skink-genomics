### Phylogeny

making a phylogenetic analysis using IQTREE.

First, I convert the SNPs to a phylogenetic input format for IQtree: phylip
by cloning the github into my directory [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip)


```
git clone https://github.com/edgardomortiz/vcf2phylip.git
cd vcf2phylip
```
since my vcf file is located just the directory above i can run ../ after the input command to run the code.

```
python vcf2phylip.py -i ../robust_mindp4_r06.recode.vcf
```
this produces robust_mindp4_r06.recode.min4.phy which will be used to submit my job 
 
I run my IQtree job requesting 24hrs, using 64G, the GTR+G model with 1000 bootstraps.

mkdir -p iqtree_output

script for job

#!/bin/bash
#SBATCH --job-name=iqtree
#SBATCH --output=iqtree_output.log
#SBATCH --error=iqtree_error.log
#SBATCH --time=24:00:00            # Max run time (hh:mm:ss)
#SBATCH --mem=64G                  # Memory per node
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1

module load IQ-TREE

iqtree2 -s robust_mindp4_r06.recode.min4.phy -m GTR+G -bb 1000 -nt 16 -pre robust_tree


```sh
module load IQ-TREE

iqtree2 -nt 16 -s robust_mindp4_r06.recode.min4.phy -st DNA -m GTR+G -bb 1000  -pre inferred
```

maybe i do iqtree -s robust_mindp4_r06.recode.min4.phy -m MFP -bb 1000 -nt AUTO which will allow me to determine best model

download FIGtree to visualise it.

# make own conda environment 
```
module load miniconda3
conda create -n iqtree_env -c bioconda iqtree
conda activate iqtree_env
```

GGtree, much much prettier (tutorial online)
