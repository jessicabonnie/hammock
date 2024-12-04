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
SCRIPT_DIR=$(dirname "$0")

if [ $# -eq 0 ]; then
    echo "Usage: $0 <input_bed> <output_prefix> [--balance] [--reuse]"
    exit 1
fi

balance_string=""
reuse=false

if [ $# -gt 2 ]; then
    if [ "$3" == "--balance" ] || [ "$3" == "--reuse" ]; then
        if [ "$3" == "--balance" ]; then
            balance_string="--balance"
        fi
        if [ "$3" == "--reuse" ]; then
            reuse=true
        fi
        
        if [ $# -gt 3 ]; then
            if [ "$4" == "--reuse" ]; then
                reuse=true
            else
                echo "Invalid option: $4"
                exit 1
            fi
        fi
    else
        echo "Invalid option: $3"
        exit 1
    fi
fi

# Check if reuse flag is not set (empty string)
# If reuse=true was not passed as an argument, execute the code below
if [ "$reuse" = false ]; then
# Generate a list of BED files of different ratio of overlap with the original
python3 $SCRIPT_DIR/generate_beds.py $inbed 2,3,10 $outname > ${outname}.list

grep -v modeB $outname.list > ${outname}A.list
grep -v modeA $outname.list > ${outname}B.list
realpath $inbed > ${outname}_primary.list
fi

parallel --jobs 4 srun -n $SLURM_CPUS_ON_NODE python3 $CODEDIR/bed_similarity.py ${outname}.list  ${outname}_primary.list --mode {1} --{2} -o ${outname} $balance_string :::  A B ::: hyperloglog minhash exact 

parallel --jobs 4 srun -n $SLURM_CPUS_ON_NODE python3 $CODEDIR/bed_similarity.py ${outname}.list  ${outname}_primary.list --mode C --{1} -o ${outname} --subsample {2} $balance_string ::: hyperloglog minhash exact :::  .1 .25 .5 .75 1



# srun python3 $CODEDIR/bed_jaccards_parallel_1m.py ${outname}.list  ${outname}_primary.list C ${outname}_real 

# parallel python3 $CODEDIR/bed_jaccards_parallel_1m.py ${outname}{1}.list  ${outname}_primary.list {1} real_${outname} ${outname} ::: A B 

