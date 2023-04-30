%% This is the main execution file.
clc; clear; close all;

%% Load the subject file to get the timezone -- This is to run data on WDAP
fid = fopen('real_data/subject.txt');
tline = fgetl(fid);
sub_file = cell(0,1);
while ischar(tline)
    sub_file{end+1,1} = tline;
    tline = fgetl(fid);
end
fclose(fid);
timezone = sub_file{2};

ref_date = datenum(1970,1,1);
bin_size = 5;

try
    if(isempty(timezone))
	t_offset=0;
    else
    	t_offset = tzoffset(datetime('today', 'TimeZone', timezone));
    	t_offset = datenum(t_offset);
    end
catch
    t_offset = 0;
end

try
    Thr = readtable('real_data/heart_rate.csv','Delimiter', ',');
    Thr = table2array(Thr(:,1:2));
    Thr(:,1) = Thr(:,1)/(24*60*60) + ref_date + t_offset;
catch
    error('No heart rate file in directory.')
end

try
    Tsteps = readtable('real_data/steps.csv', 'Delimiter', ',');
    Tsteps = table2array(Tsteps(:,1:2));
    Tsteps(:,1) = Tsteps(:,1)/(24*60*60) + ref_date + t_offset;
    isEmptySteps = 0;
catch
    Tsteps = [];
    isEmptySteps = 1;
end

try
    Tsleep = readtable('real_data/sleep.csv', 'Delimiter', ',');
    Tsleep = table2array(Tsleep(:,1:2));
    Tsleep(:,1) = Tsleep(:,1)/(24*60*60) + ref_date + t_offset;
    isEmptySleep = 0;
catch
    Tsleep = [];
    isEmptySleep = 1;
end 

if(~isEmptySleep)
    Thr = remove_sleep(Thr, Tsleep);
end
if(~isEmptySteps)
    %% Run Heart rate Algorithm (Bowman et al., 2021)
    dmy_hr = bayes_hr_estimator(Thr, Tsteps, bin_size);
    
    %% Execute the main algorithm
    LSKF_circadian(dmy_hr, Thr, Tsteps, bin_size);
end