##############
# this shell scripts allow to loop over the inputs CSV file in : inp_path
##############

# inputs files directory and list: mind!! that we sue the structure of the input-file-name to give
# correct argval to python command ---> see below 

inp_path=~/aerosol/timo_input/
ii_files=`ls -v $inp_path*i_to_i*`
ir_files=`ls -v $inp_path*i_to_r*`

 
arr_ii=($ii_files)
arr_ir=($ir_files)

num_cond=${#arr_ir[@]}

# echo ${#arr_ii[@]}
 

# parameters of the simulation
#INDEX_CASE=1
#RUNS=100

# loop over inputs files
 
#for i in  $(seq 15 $num_cond); do # here we forgot to add the 0th element!!!!!: (seq 0 $num_cond)  
for i in  $(seq 14 14); do
 		
	f_ii=${arr_ii[i]}
	f_ir=${arr_ir[i]}

	# get the name of condition from the input file

	filename="${f_ii##*/}"  
	OUTPUT_SIGN1=${filename#i_to_i_}
	OUTPUT_SIGN2=${OUTPUT_SIGN1%_schedule.csv}
	echo $OUTPUT_SIGN2

	if [ ! -d "$OUTPUT_SIGN2" ]; then   # check if dir. exists
		if [ ! -L "$OUTPUT_SIGN2" ]; then # if It is a symlink!...
	  		mkdir $OUTPUT_SIGN2	 
  	  	fi
	fi
											 
	cd $OUTPUT_SIGN2
	
	time python2.7 ~/aerosol/aerolsol_model_outputs/aerosol_model_cpu_par.py \	
	~/aerosol/timo_input/school_i_j_aggreg.txt \
	$f_ii $f_ir $OUTPUT_SIGN2

	

	cd ../
done
