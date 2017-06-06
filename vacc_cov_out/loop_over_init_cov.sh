# correct argval to python command ---> see below 
 
f_ii=~/aerosol/timo_input/i_to_i_ACR05_schedule.csv
f_ir=~/aerosol/timo_input/i_to_r_ACR05_schedule.csv

# get the name of condition from the input file
filename="${f_ii##*/}"  
OUTPUT_SIGN1=${filename#i_to_i_}
OUTPUT_SIGN2=${OUTPUT_SIGN1%_schedule.csv}
printf 'school condition: %s \n' "$OUTPUT_SIGN2"


# loop over alpha values : NB!: this interger numb alpha_fac is devided by
# 10 in the .py script ----> ALPHA = alpha_fac/10.0

for coverage_fract in  $(seq 0 10 100); do 


	if [ ! -d "$coverage_fract" ]; then   # check if dir. exists
		if [ ! -L "$coverage_fract" ]; then # if It is a symlink!...
  	  		mkdir $coverage_fract	 
  	  	fi
	fi
	
	cd $coverage_fract

 	printf 'computing sim. for vaccination coverage = %d %%\n' "$coverage_fract"

	time python2.7 ~/aerosol/vacc_out/vac_fully_aerosol_cpu_par.py \
	$f_ii $f_ir $OUTPUT_SIGN2 $coverage_fract
	
	cd ../
done
