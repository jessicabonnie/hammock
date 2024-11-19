#! /bin/bash
#SBATCH --job-name=parallel_srun
#SBATCH --output=parallel_srun_%j.out
#SBATCH --error=parallel_srun_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G  # Increased memory allocation
#SBATCH --partition=parallel

ml parallel
ml anaconda
conda activate minhash
export PYTHONPATH=$PYTHONPATH:/home/jbonnie1/interval_sketch/

inbed=$1 #../../remap2022_hg38_1M.bed
outname=$2 #remap1M
CODEDIR=/home/jbonnie1/interval_sketch/hammock/lib
SLURM_CPUS_ON_NODE=${SLURM_CPUS_ON_NODE:-1}

python3 $CODEDIR/generate_beds.py $inbed 2,3,4,5 $outname > ${outname}.list

grep -v modeB $outname.list > ${outname}A.list
grep -v modeA $outname.list > ${outname}B.list
realpath $inbed > ${outname}_primary.list


parallel --jobs 4 srun -n $SLURM_CPUS_ON_NODE python3 $CODEDIR/bed_similarity.py ${outname}.list  ${outname}_primary.list --mode C --{1} -o ${outname} --subsample {2} ::: hyperloglog minhash exact :::  .1 .25 .5 .75 1
# for subsample in .1 .25 .5 .75 1; do
# echo "Subsample: $subsample"
#     for perm in 50 100 200 500 ; do
#         echo "Permutation number: $perm"
#         srun -n $SLURM_CPUS_ON_NODE python3 $CODEDIR/bed_jaccmh_parallel_1m.py ${outname}.list  ${outname}_primary.list --mode C --perm $perm -o ${outname} --subsample $subsample
# echo "Complete: Permutation number: $perm "
#     done
# done



# parallel --jobs 4 srun -n $SLURM_CPUS_ON_NODE python3 $CODEDIR/bed_jaccmh_parallel_1m.py ${outname}{1}.list  ${outname}_primary.list --mode {1} --perm {2} -o ${outname} ::: A B ::: 50 100 200 500

# srun python3 $CODEDIR/bed_jaccards_parallel_1m.py ${outname}.list  ${outname}_primary.list C ${outname}_real 

# parallel python3 $CODEDIR/bed_jaccards_parallel_1m.py ${outname}{1}.list  ${outname}_primary.list {1} real_${outname} ${outname} ::: A B 


# for perm in 50 100 200 500 ; do
# for mode in "A" "B"; do
#  echo "Processing mode: $mode"
#     python3 $CODEDIR/bed_jaccmh_parallel_1m.py ${outname}${mode}.list  ${outname}_primary.list --mode $mode --perm $perm -o ${outname}
#     done
# done
