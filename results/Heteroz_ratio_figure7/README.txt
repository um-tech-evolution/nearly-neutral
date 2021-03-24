CSV and jl files copied from:  evotech/nearly-neutral/plots/9_13_17fixed.

The file: in_fixed_mulist_N25_6400_s_all_ratios.csv is used to create the plot.

This file was created "by hand" by concatenating the CSV files for the various s values.

3/17/21:  What are the raios and How are they computed?  
Ratios are a selection coefficient quantity divided the corresponding neutral quantity.
For the 2017 paper the ratios were computed in the spreadsheet: in_fixed_mulist_N25_6400_s_all_ratios.xlsx
computed_count_ratio = ave_innov_count/neutral_ave_innov_count where neutral_ave_innov_count is the neutral run ave_innov_count
computed_obs_count_ratio = ave_heteroz/neutral_ave_heteroz
computed_exp_count_ratio = 2*(4*C28*B28-1+EXP(-4*C28*B28))/(4*C28*B28*(1-EXP(-4*C28*B28)))  for row 28
    where column B is s and column C is N
computed_exp_count_ratio = 2(4*N*s-1+exp(-4*N*s))/(4*N*s*(1-exp(-4*N*s)))

Finally, 
count_ratio = computed_count_ratio with 1 for the neutral rows (except not computed for some rows)
observed_heteroz_ratio = computed_obs_count_ratio with 1 for the neutral rows (except not computed for some rows)
expected_heteroz_ratio = copied from computed_exp_count_ratio but seems to be off by 1 row
