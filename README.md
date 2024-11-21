######################################################################################## 
SV-Quest 2.0

A Perl scripts to call indel position and probable IS transposition from mapped.bam.   

SV Quest: Sensitive Structure Variation detection tool

Kazuma Uesaka, Yuichi Fujita, and Kunio Ihara  



Input: 
  bam file and reference.fasta for the mapping   
Outnput:	
  indel and deletion position printed to STDOUT  

Usage:  
  SV-Quest_v2.2.6.pl -f reference.fa -1 forward.fq -2 reverse.fq -b IS.fasta -n sample_prefix -y 4

########################################################################################


<p align="center"><img src="Figure1.png" alt="workflow" width="800"></p>


    
## 1. Requirements  
- SAMTools  (version >= 1.3.1)  
- BWA (version >= 0.7.17)  
- sambamba  (version >= 0.8.0)
- SPAdes  (version >= 3.14)
- Minimap2
- Mirabait  



Install Dependancy using conda package manager (macos, Linux(WSL)).  

```
mamba create -n SV-Quest2 python=3.10 -y
conda activate SV-Quest2
mamba install -c bioconda bwa samtools sambamba spades blast minimap2
#mirabait (using apt package maneger or download from Mirabait repositoru)
sudo apt upgrade and apt install mira-assembler
```
    


## 2. clone source
```
cd $HOME/bin/ 
git clone https://github.com/kazumaxneo/SV-Quest2.git
cd SV-Quest2/
chmod +x  SV-Quest_v2.2.6.pl
echo export PATH=\$PATH:`pwd`\ >> ~/.bash_profile && source ~/.bash_profile
SV-Quest.pl
```
  

  
## 3. Run  
1, prepare Tn sequernces.  
I recommned isescan program to identify endogenous IS sequernces from ref. genome. 
```
mamba install -c bioconda isescan -y
isescan.py --seqfile reference.fna --output results --nthread 8
```
The output directory xxx.fna.orf.fna is the DNA sequence of the IS detected.  

2, Run SV-Quest2.  
The IS sequence is specified with -b </path/to/file>.  
The reference genome is specified with -f </path/to/file>. 
```
SV-Quest_v2.2.6.pl -f reference.fa -1 forward.fq -2 reverse.fq -b IS.fasta -n sample1 -y 4
```
  


 
## Test Run 
```
git clone https://github.com/kazumaxneo/SV-Quest2.git
cd SV-Quest2/
mkdir test && cp ISs.fna test/
#fastq and ref.fasta file
git clone https://github.com/kazumaxneo/SV-Quest.git
tar zxvf SV-Quest/sample.tar.gz -C test/
cd test/
SV-Quest_v2.2.6.pl -f chromosome.fasta -1 forward.fq -2 reverse.fq -b ISs.fna -n sample1 -y 4
```

## Licence ##

GPL v3.


