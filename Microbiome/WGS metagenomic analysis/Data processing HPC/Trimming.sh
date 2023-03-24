### This script trims low/poor quality reads from a fastq file ###
### To run this script type "sh Trimming.sh" inside the fastq directory ###
for i in /gpfs/projects/whendrickx_Lab/SDR100029/WGS/WGS_all/unaligned_fastq/*_unaligned.fastq; do 
name=$(basename "$i" _unaligned.fastq);
if [ -e "/gpfs/projects/mharis_Lab/Data/microbiome/trimmed/"${name}"_unaligned_trimmed.fq" ]
 then
  echo "fastq file exists"
 else
  echo " fastq file added"
  B1="/gpfs/projects/whendrickx_Lab/SDR100029/WGS/WGS_all/unaligned_fastq/"$name"_unaligned.fastq"
  ID="${name}"
  echo "sh Trimming_trimgalore.sh $B1 /gpfs/projects/mharis_Lab/Data/microbiome/trimmed" $ID | bsub -n 8 -e /gpfs/projects/mharis_Lab/Data/microbiome/trimmed/logs/${name}.err -o /gpfs/projects/mharis_Lab/Data/microbiome/trimmed/logs/${name}.out -P TRIM
fi
done