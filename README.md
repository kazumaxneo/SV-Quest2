######################################################################################## 
SV-Quest 2.0

A Perl scripts to call indel position from mapped.bam.   

SV Quest: Sensitive Structure Variation detection tool

Kazuma Uesaka, Yuichi Fujita, and Kunio Ihara  



Input: 
  bam file and reference.fasta for the mapping   
Outnput:	
  indel and deletion position printed to STDOUT  

Usage:  
  perl SV-Quest.pl -f reference.fa -1 forward.fq -2 reverse.fq -b transposon.fasta -n sample1 -y 4


 The mapped.bam and it's reference.fasta should be included in the temp folder,  
 and copy the Indel_Hunter.pl in this folder.
########################################################################################


<p align="center"><img src="Figure1.png" alt="workflow" width="800"></p>


    
## Requirements  
- SAMTools  (version >= 1.3.1)  
- BWA (version >= 0.7.17)  
- sambamba  (version >= 0.8.0)
- SPAdes  (version >= 3.14)  



Install Anaconda (Mac OS X, Linux).  

```
mamba create -n SV-Quest2 python=3.10 -y
conda activate SV-Quest2
mamba install -c bioconda bwa samtools sambamba spades

```
    


## Source
```
cd $HOME/bin/ 
git clone https://github.com/kazumaxneo/SV-Quest2.git
cd SV-Quest2/
chmod +x  SV-Quest_v2.2.6.pl
echo export PATH=\$PATH:`pwd`\ >> ~/.bash_profile && source ~/.bash_profile
SV-Quest.pl
```
    


## Licence ##

GPL v3.


