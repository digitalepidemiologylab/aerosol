# correct argval to python command ---> see below 
 
f_ii=/home/lazzari/school_net_proj/input_files/i_to_i_ACR05_schedule.csv
f_ir=/home/lazzari/school_net_proj/input_files/i_to_r_ACR05_schedule.csv
# get the name of condition from the input file
filename="${f_ii##*/}"  
OUTPUT_SIGN1=${filename#i_to_i_}
OUTPUT_SIGN2=${OUTPUT_SIGN1%_schedule.csv}
printf 'school condition: %s \n' "$OUTPUT_SIGN2"



for coverage_fract in  $(seq 0 10 100); do 

	if [ ! -d "$coverage_fract" ]; then   # check if dir. exists
		if [ ! -L "$coverage_fract" ]; then # if It is a symlink!...
  	  		mkdir $coverage_fract	 
  	  	fi
	fi
	
	cd $coverage_fract

 	printf 'computing sim. for vaccination coverage = %d %%\n' "$coverage_fract"

	time python2.7 ~/aerosol/1000_comb+vacc/vacc+combined_model_cpu_par.py \
	~/aerosol/input_files/school_i_j_aggreg.txt \
	$f_ii $f_ir $coverage_fract $OUTPUT_SIGN2
	
	cd ../
done
