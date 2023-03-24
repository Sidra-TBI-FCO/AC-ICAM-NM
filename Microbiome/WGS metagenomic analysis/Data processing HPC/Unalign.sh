### This script extracts unaliogned reads from a WGS/WES bam file ###
### To run this script type "sh Unalign.sh" inside the bam directory ###
for i in BAMS/*.sorted.bam; do 
name=$(basename "$i" .sorted.bam);
if [ -e "unaligned_bams/"${name}".unaligned.sorted.bam" ]
 then
  echo "bam file exists"
 else
  echo " bam file added"
  B1="BAMS/"$name".sorted.bam"
  ID="${name}"
  echo "sh Unalign_Samtools.sh $B1 /gpfs/projects/whendrickx_Lab/SDR100029/WGS/WGS_all/unaligned_bams" $ID | bsub -n 8 -e unaligned_bams/logs/${name}.err -o unaligned_bams/logs/${name}.out -P UNALIGN
fi
done