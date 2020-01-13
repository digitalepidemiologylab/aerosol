
Code and data for the paper 
["Assessing the Dynamics and Control of Droplet- and Aerosol-Transmitted Influenza Using an Indoor Positioning System"](https://www.nature.com/articles/s41598-019-38825-y)


- folder [location_info](https://github.com/salathegroup/aerosol/tree/master/location_info) contains the original data from the experiment ([Salathé et al., 2010](http://www.pnas.org/content/107/51/22020.short)) and the code used to generate the simulated exposures 

- folder [input_files](https://github.com/salathegroup/aerosol/tree/master/input_files) contains the inputs needed to run the simulations: one text (.txt) file for the recorded CPIs and table (.csv) files for the simulated quanta and exposures, for each ventilation scenario

- folders [vacc_cov_out](https://github.com/salathegroup/aerosol/tree/master/vacc_cov_out) and [aerolsol_model_outputs](https://github.com/salathegroup/aerosol/tree/master/aerolsol_model_outputs) contain the code to simulate outbreak in the _full aerosol_ scenario, with and without vaccination

- folders [1000_comb](https://github.com/salathegroup/aerosol/tree/master/1000_comb) and [1000_comb+vacc](https://github.com/salathegroup/aerosol/tree/master/1000_comb%2Bvacc) contain respectively the code to simulate outbreak in the _combined-model_ scenario, with and without vaccination

- folder [analysis_simul_out](https://github.com/salathegroup/aerosol/tree/master/analysis_simul_out) contains the Rmd notebooks (also in html format) used to produce some figures useful for analysis, including fig. 3 and fig. 4 of the paper

- folder [r0](https://github.com/salathegroup/aerosol/tree/master/r0) contains the python code (.py), the files the run the simulation (.sh) and results (.csv) used to compute R'0, for the three different models (aerosol, droplet, combined)

- folder [classroom_fig](https://github.com/salathegroup/aerosol/tree/master/classroom_fig) contains the data and the code to reproduce fig. 2 of the paper

