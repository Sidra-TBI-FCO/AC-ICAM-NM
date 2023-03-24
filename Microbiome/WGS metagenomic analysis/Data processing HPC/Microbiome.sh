### This script performs metagenomics analysis of fastq file ###
### and outputs microbiome relative abundances at different bacterial levels  ###
### the script can be modified to run for other kindoms like viruses, archea etc by adjusting --ignore option ###
### To run this script type "sh Microbiome.sh" inside the fastq directory ###

for i in trimmed/*_trimmed.fq; do 
name=$(basename "$i" _trimmed.fq);
if [ -e "aligned/"${name}"_a.txt" ]
 then
  echo "file exists"
 else
  echo " file added"
  B1="trimmed/"$name"_trimmed.fq"
  ID="${name}"
  echo "sh Microbiome_metaphlan.sh $B1 /gpfs/projects/mharis_Lab/Data/microbiome/aligned" $ID | bsub -n 8 -e /gpfs/projects/mharis_Lab/Data/microbiome/aligned/logs/${name}.err -o /gpfs/projects/mharis_Lab/Data/microbiome/aligned/logs/${name}.out -P MICROBIOME
fi
done