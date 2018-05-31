# inputs files directory
# inp_path=/home/lazzari/school_net_proj/input_files/*.csv
 

time python2.7 ../r0_vacc+droplet_alpha0_model_cpu_par.py \
../../input_files/school_i_j_aggreg.txt \
../../input_files/i_to_i_ACR05_schedule.csv \
../../input_files/i_to_r_ACR05_schedule.csv \
&&  echo  finished at `date`
