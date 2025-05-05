# *Sexing 
load up miniconda 
```
module load Miniconda3/23.10.0-1
```

conda plays up so have to purge and reload in a specific way

```
module purge && module load Miniconda3
```

```
source $(conda info --base)/etc/profile.d/conda.sh
```
```
export PYTHONNOUSERSITE=1
```

might not need conda config part
```
conda config --add pkgs_dirs /nesi/nobackup/uoo04250/$USER/conda_pkgs
```
```
conda create -n radsexn
```

```
conda activate radsexn
```

```
radsex
```

```
radsex process --input-dir ./samples --output-file markers_table.tsv --threads 16 --min-depth 1
radsex distrib --markers-table markers_table.tsv --output-file distribution.tsv --popmap popmapbaldefg.tsv --min-depth 5 --groups M,F
```


```
#nbootstraps=500
file_sample="popmapbaldefg.tsv"
outfile="500radsexbootstrap.txt"
#radsex process --input-dir ./samples --output-file markers_table.tsv --threads 10 --min-depth 1
rm $outfile
for i in {1..500}
do
    echo $i
    #shuffle table
    cut -f 1 $file_sample | shuf | paste -d "\t"  - $file_sample  | cut -f 1,3 > tempsex.txt
    radsex distrib --markers-table markers_table.tsv --output-file tempdistrib.txt     --popmap tempsex.txt --min-depth 1 --groups M,F
     grep -E "^2\s+15"  tempdistrib.txt | cut -f 3 > temponlyf.txt
 grep -E "^0\s+16"  tempdistrib.txt  | cut -f 3 >> temponlyf.txt
  grep -E "^1\s+16\s+"  tempdistrib.txt  | cut -f 3 >> temponlyf.txt
  grep -E "^2\s+16"  tempdistrib.txt  | cut -f 3 >> temponlyf.txt
  grep -E "^0\s+15\s+"  tempdistrib.txt  | cut -f 3 >> temponlyf.txt
  grep -E "^1\s+15\s+"  tempdistrib.txt  | cut -f 3 >> temponlyf.txt
awk '{sum += $1} END {print sum}' temponlyf.txt > femalesum.txt
     grep -E "^15\s+2"  tempdistrib.txt  | cut -f 3 > temponlym.txt
  grep -E "^16\s+0"  tempdistrib.txt  | cut -f 3 >> temponlym.txt
  grep -E "^16\s+1\s+"  tempdistrib.txt  | cut -f 3 >> temponlym.txt
  grep -E "^16\s+2"  tempdistrib.txt  | cut -f 3 >> temponlym.txt
  grep -E "^15\s+0\s+"  tempdistrib.txt  | cut -f 3 >> temponlym.txt
  grep -E "^15\s+1\s+"  tempdistrib.txt  | cut -f 3 >> temponlym.txt
awk '{sum += $1} END {print sum}' temponlym.txt > malesum.txt
    paste femalesum.txt malesum.txt >> $outfile
    #rm tempsex.txt temponlyf.txt temponlym.txt tempdistrib.txt # Remove table
done
```
allows to test combination markers with the bootstrap analysis focused on the most extreme sex-biased regions of the distribution (combinations with 0-2 individuals of one sex and 15-16 of the opposite sex) as these regions are the most informative when inferring sex-linked markers.
