# Level set Kalman Filter for estimation of the circadian phase and its uncertainty from wearable data
### Developers: Dae Wook Kim, Minki P. Lee, Daniel B. Forger
### Overview: This algorithm was developed to analyze the underlying circadian rhythm in the master circadian clock using wearable heart rate, steps, and sleep data.

## Attached file lists
### Computer codes
main.m: This is the main execution file for analyzing the wearable heart rate and activity data to estimate the circadian phase and its uncertainty of the master circadian pacemaker.

bayes_hr_estimator.m: This function applies the Monte-Carlo Markov chain method to wearable data to obtain the circadian phase and uncertainty of heart rate circadian clock (Bowman et al., 2021). 

LSKF_circadian.m: This is the main function that does preprocessing, sets up for the implementation of the Level Set Kalman Filter (LSKF), and saves the result.

lskf.m: This is the function that implements the LSKF.

avg_vel_time_update.m: This is the function that implements forward propagation of the circadian clock system and estimate the predicted circadian phase, using the level set method (Wang and Forger, 2021).

simpler_circadian.m: This is the function that definees the velocity field of the system based on the mathematical model of the master circadian clock (Forger et al., 1999).

measurement_update.m: This is the function that implements the measurement-update step based on the square root cubature method (Arasaratnam et al., 2009).

### Wearable data
subject.txt: A basic subject information including deidentified subject ID and time-zone.

heart_rate.csv: A two column array, the first column lists the epoch time date in days (e.g., 738157 is 01-Jan-2021), and the second colum lists the heart rate value at that time.

steps.csv: A two column array like the heart rate data, but the second column is the steps value at the respective time.

sleep.csv: A two column array like the heart rate data, but the second column describes whether the subject sleep. Zero value denotes wakefulness, and positive values denote sleep. Even if we have no sleep data, we can use the algorithm only with heart rate and steps data.

### Output files
visualiation: Visualization of the level set ellipsoids on x,xc plane are plotted for each day. Blue ellipsoids correspond to the propagation before the correction (i.e., during the time-update step). Red ellipsoids correspond to the propagation after the correction (i.e., after the measurement-update step). Set varialbe "plotLevelSet" as 0 to turn off visaulizations.

hr_phase_result.csv: This is the .csv output file from LSM_hr_estimator_v9.m that contains information about heart rate phase, heart rate phase uncertainty, and heart rate circadian rhythm paramters.

result.csv: This is the main result file that from the LSKF_circadian.m. Each column corresponds to the date of the data collection. The columns contain information about the 
• corrected circadian phase and uncertainty of the master clock (i.e., after the measurement-update step) 
• predicted circadian phase and uncertainty of the master clock (i.e., before the measurement-update step) 
• phase estimates and uncertainty of the heart rate circadian clock (same as in hr_phase_result.csv)
• estimated circadian phase and uncertainty of the master clock solely based on the mathematical model and level set method (i.e., result without measurement-update step).
