
choose=$1
node=$2
node_ppn=1
from=$3
to=$4
for sb in $(seq $from $to);do
	chr=$sb
	choose=$choose
	process=$choose"_"$chr
	binsize=200

	root_path="/share/data6/tmp/renjun/ExtraJobs/202206_HuR_RIP-Seq"
	script_path=$root_path/1_mainSeek.R
	run_job=$root_path/7_RIPSeeker/$process.job

	### Rscript job --- 59
	# echo "#! /bin/bash" > $run_job
	# echo "#SBATCH --job-name=$process" >> $run_job
	# echo "#SBATCH -w $node" >> $run_job
	# echo "#SBATCH --time 20000:00:00" >> $run_job
	# echo "#SBATCH --nodes 1" >> $run_job
	# echo "#SBATCH --ntasks 1" >> $run_job
	# echo "#SBATCH --output=$process_%j.log" >> $run_job
	# echo "source activate R4.1 " >> $run_job
	# echo "Rscript $script_path $chr $choose $binsize" >> $run_job
	# sbatch $run_job

	### Rscript job --- 121
	echo "#! /bin/bash" > $run_job
	echo "#PBS -N $process" >> $run_job
	echo "#PBS -l nodes=node$node:ppn=$node_ppn,walltime=20000:00:00" >> $run_job
	echo "#PBS -j oe" >> $run_job
	echo "source activate R4.1 " >> $run_job
	echo "Rscript $script_path $chr $choose $binsize $root_path" >> $run_job
	qsub $run_job

done
