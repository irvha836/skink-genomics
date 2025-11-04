first need all the dependancies for snailpot generation by cloning this GitHub into my directory 
https://github.com/rjchallis/assembly-stats

```
module load Perl

perl pl/asm2stats.minmaxgc.pl /nesi/nobackup/uoo04250/genome_assembly/psmcclean/whitakermasked/whitakergenome.fasta > whitaker_assembly2.minmaxgc.json
```
tagged on repeat content and busco onto 


downloaded all the files from the repository onto my own computer 
updated the HTML file to point to my own json file 

ran this to override chrome loading issues
```
python3 -m http.server 8000
```
plot with this 

```
http://localhost:8000/assembly-stats.html
```
