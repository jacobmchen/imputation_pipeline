# this file will preprocess the biomarker data to create a new column
# representing the average of each biomarker's change compared to baseline
# at the 4 time measured time points of each biomarker
# in math, (time_1 - time_0 + time_2 - time_0 + time_3 - time_0 + time_4 - time_0)/4
# = (time_1 + time_2 + time_3 + time_4 - 4*time_0)/4
# = (time_1+time_2+time_3+time_4)/4 - time_0
