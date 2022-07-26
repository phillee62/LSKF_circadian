function dmy_hr = LSM_hr_estimator_v9(hr_data, steps_data, bin_size)

%% Close figures

close all

%% Inputs:

% 1. hr_data - a two column array, the first column lists the epoch
% time date in days (e.g., 738157 is 01-Jan-2021), the second column lists
% the heart rate value at that time
% 2. steps_data - a two column array like the heart rate data, but the
% second column is the steps value at the respective time
% 3. bin_size - the number of minutes to bin the data (5 is suggested)

%% Outputs:

% A results.csv file that lists the parameter estimates for each day there
% is data along with the uncertainties in the estimates and the number of
% data points for that day and the daily step count for that day

%% Global parameters
global hr_avg step_avg N bin_size1

bin_size1 = bin_size;

%% Initialize parameters

% Credible interval should contain this fraction of samples
% confidence = 0.8;
% nbars = 48; % Number of bars per day (so 48 is half-hour activity bars)
% thresh = 250; % Number of steps in bin for bar to be full
% minbar = 30; % Bins with fewer than this number of steps aren't shown
% nytick = 15; % Tick every this number of days on the Y Axis

% Optimal parameter values and uncertainties are stored in these lists

%all_fits = [];
% all_stds = [];
%phases = [];
%HR_phase_sample_collection = [];

%% Open the results file

% Open results.csv file and add a header line
results_file = fopen('hr_phase_results.csv', 'w');
% results_file = fopen('results.csv', 'w');
fprintf(results_file, 'Date (dd-mmm-yyyy), Horizontal Shift (Phase Estimate-Raw), Phase Estimate (Time), STD of Phase Estimate (hr), Number of data points, Daily step count \n');
%% Process data
hr_dates = hr_data(:,1);
steps_dates = steps_data(:,1);

hr_vals = hr_data(:,2);
steps_vals = steps_data(:,2);

% Save all the days in the interval of data given to loop over
hr_days = floor(hr_dates);
dmy_hr = unique(hr_days, 'rows', 'stable');
dmy_hr = (dmy_hr(1):(dmy_hr(end)+1));
num_days = length(dmy_hr)-2;

%% Remove hr_times and hr_values where steps are zero consecutively for more than 4 hours <- 4 hours might be typo... 2 hours would be correct
[hr_dates, hr_vals] = constant_steps_filter2(hr_dates, hr_vals, steps_dates, steps_vals,120);

%% Initialize the figure
% num_mins = 24/nbars*60;
% f = figure(1);
% xlim([0 48])
% ylim([1 num_days+1])
% set(gca, 'YDir','reverse')
% hold on

%% Main loop

epara=0;

for i = 1:num_days
    
    % Load steps information for that day
    steps_times = steps_dates(steps_dates>=dmy_hr(i) & steps_dates < dmy_hr(i+2));
    steps_values = steps_vals(steps_dates>=dmy_hr(i) & steps_dates < dmy_hr(i+2));
    
    % Load heart rate data for the specific day
    hr_times = hr_dates(hr_dates < dmy_hr(i+2) & hr_dates >= dmy_hr(i));
    hr_values = hr_vals(hr_dates < dmy_hr(i+2) & hr_dates >= dmy_hr(i));
    
    % Convert the times data to minutes
    hr_times = hr_times*24*60;
    steps_times = steps_times*24*60;
    
    daily_step = sum(steps_values);
    
    % Ask whether the hr data is empty for that specific day. If so, make
    % the fit values NaN and continue in the loop
    
    if(isempty(hr_times))
        
        %all_fits = [all_fits; NaN];
        %all_stds = [all_stds; 0];
        %if(i==1)
        %    phases =[];
        %else
        %    phases =[phases; phase_est];
        %end
        
        fprintf(results_file, '%s, %f, %s, %f, %i, %i\n', datestr(dmy_hr(i+1)), NaN, NaN, NaN, NaN, NaN);
        continue;
        
        
    else
        
        raw_steps_data = [steps_times steps_values];
        raw_hr_data = [hr_times hr_values];
        
        % Find all unique binned times for HR data
        [hr_avg, ~, idx] = unique(floor(raw_hr_data(:, 1) / bin_size), 'stable');
        
        % Average HR measurements in same bin
        val = accumarray(idx, raw_hr_data(:, 2), [], @mean);
        
        % Combine unique times and averages
        hr_avg = [hr_avg, val];
        
        try
            
            % Find leftmost and rightmost bin
            left_min = floor(min(raw_hr_data(1, 1), raw_steps_data(1, 1)) / bin_size);
            right_max = floor(max(raw_hr_data(end, 1), raw_steps_data(end, 1)) / bin_size);
            clearvars raw_hr_data;
            
        catch
            
            left_min = floor(raw_hr_data(1, 1) / bin_size);
            right_max = floor(raw_hr_data(end, 1) / bin_size);
            clearvars raw_hr_data;
            
        end
        
        % Total length of interval containing HR data
        period_offset = [left_min, right_max - left_min + 1] * bin_size;
        
        % Precomputed quantities for efficiency: gaps in consecutive measurements
        %jumps = [hr_avg(1, 1) + 0.1; hr_avg(2:end, 1) - hr_avg(1:(end - 1), 1)];
        
        % " " " ": number of measurements after averaging same times
        N = size(hr_avg, 1);
        
        % Fill in any gaps in step data with zeros
        step_int = int32(raw_steps_data(:, 1));
        steps_new = zeros(period_offset(2), 2);
        steps_new(:, 1) = period_offset(1) + (0:(period_offset(2) - 1));
        steps_new(step_int - period_offset(1) + 1, 2) = raw_steps_data(:, 2);
        clearvars raw_steps_data;
        
        % Average steps data into bins (including zeros)
        [step_avg, ~, idx] = unique(floor(steps_new(:, 1) / bin_size) , 'stable');
        val = accumarray(idx, steps_new(:, 2), [], @mean);
        step_avg = [step_avg, val];
        
        % We only need to keep the step bins corresponding to HR data
        step_avg = step_avg(hr_avg(:, 1) - step_avg(1, 1) + 1, :);
        
        if(length(hr_avg(:,1))<20)
            
            %all_fits = [all_fits; NaN];
            %all_stds = [all_stds; 0];
            %phases = [phases; phase_est];
            fprintf(results_file, '%s, %f, %s, %f, %i, %i\n', datestr(dmy_hr(i+1)), NaN, NaN, NaN, NaN, NaN);
            continue;
            
        end
        
        time_diff=diff(hr_avg(:,1));
        data_dim1 = length(time_diff);
        heart_rate1 = -100*ones(data_dim1,1); heart_rate2 = -100*ones(data_dim1,1);
        activity1 = -100*ones(data_dim1,1); activity2 = -100*ones(data_dim1,1);
        time1 = -100*ones(data_dim1,1); time2 = -100*ones(data_dim1,1);
        
        for ii = 1:data_dim1
            
            if abs(time_diff(ii)-1) < 0.01
                heart_rate1(ii) = hr_avg(ii,2);
                heart_rate2(ii) = hr_avg(ii+1,2);
                activity1(ii) = step_avg(ii,2);
                activity2(ii) = step_avg(ii+1,2);
                time1(ii) = hr_avg(ii,1) * bin_size;
                time2(ii) = hr_avg(ii+1,1) * bin_size;
                
            end
            
        end
        
        heart_rate1(heart_rate1==-100)=[]; heart_rate2(heart_rate2==-100)=[];
        activity1(activity1==-100)=[]; activity2(activity2==-100)=[];
        time1(time1==-100)=[]; time2(time2==-100)=[];
        
        data_dim2 = length(heart_rate2);
        
        if data_dim2 < 20
            fprintf(results_file, '%s, %f, %s, %f, %i, %i\n', datestr(dmy_hr(i+1)), NaN, NaN, NaN, NaN, NaN);
            continue
            
        else
            
            tau = 24; % free-running period
            
            xmatrix = zeros(data_dim2,5);
            xmatrix(:,1) = cos((2*pi/tau)*mod((time2 / 60), 24));
            xmatrix(:,2) = sin((2*pi/tau)*mod((time2 / 60), 24));
            xmatrix(:,3) = heart_rate1;
            xmatrix(:,4) = activity2;
            xmatrix(:,5) = -activity1;
            
            % creation of initial guess for nonlinear regression
            if epara == 0
                
                %ixmatrix = [ones(data_dim2,1), xmatrix];
                %initial_guess = (((ixmatrix')*(ixmatrix))\(ixmatrix'))*heart_rate2;
                %initial_guess = initial_guess(1:5)';
                initial_guess = [73, 4*cos((5/4)*pi), 4*sin((5/4)*pi), 0.8, 0.3];
                
            else
                
                initial_guess = epara';
                
            end
            
            % nonlinear regression model
            %b(1): mean heart rate * (1- correlation noise coefficient)
            %b(2): coefficient of first-order cosin term * (1- correlation noise coefficient)
            %b(3): coefficient of first-order sin term * (1- correlation noise coefficient)
            %b(4): correlation noise coefficient
            %b(5): heart rate increases per one step
            modelfun = @(b,x)b(1) + b(2)*x(:,1) + b(3)*x(:,2) + ...
                b(4)*x(:,3) + b(5)*x(:,4) +b(4)*b(5)*x(:,5);
            
            % nonlinear regression
            estimate = fitnlm(xmatrix, heart_rate2, modelfun, initial_guess);
            
            % estimation of circadian parameters and noise parameters
            eparainfo = estimate.Coefficients;
            eparainfo = eparainfo{:,:};
            
            epara = eparainfo(1:5,1)';
            cov_matrix=estimate.CoefficientCovariance;
            
            % Parameter sampling
            HR_sample_num =100000;
            sample_num=HR_sample_num; % the number of parameter samples (i.e., iteration number of MCMC)
            samples=mvnrnd(epara, cov_matrix, sample_num); % Here, normal assumption is used.
            
            parameter1 = samples(:,4); % correlation parameter
            parameter2 = samples(:,5); % heart rate increase per step
            
            % If the sampled parameter value is less than 0, it is just
            % set to 0.
            parameter1(parameter1 < 0) = 0;
            parameter2(parameter2 < 0) = 0;
            
            samples(:,4) = parameter1;
            samples(:,5) = parameter2;
            
            tsamples=[samples(:,1)./(1-samples(:,4)), samples(:,2)./(1-samples(:,4)), ...
                samples(:,3)./(1-samples(:,4)), samples(:,4), samples(:,5)];
            
            % amplitude samples
            amp_samples = sqrt(sum(tsamples(:,2:3).^2, 2));
            
            % phase samples
            cosvalues = tsamples(:,2) ./ amp_samples;
            sinvalues = tsamples(:,3) ./ amp_samples;
            phase_samples = zeros(length(cosvalues),1);
            
            for kk = 1:length(cosvalues)
                
                sol = angleCalc(sinvalues(kk), cosvalues(kk),'rad');
                sol = sol + pi;
                
                while sol < 0
                    
                    sol = sol + 2*pi;
                    
                end
                
                while sol >= 2*pi
                    
                    sol = sol - 2*pi;
                    
                end
                
                phase_samples(kk) = sol;
                
            end
            
            phase_mean = angle(sum(exp(1i*phase_samples)));
            phase_mean = mod(phase_mean, 2*pi);
            phase_std = circ_std(phase_samples);
            
            %phase_sample_store = (24/(2*pi)) * mod(phase_samples,2*pi);
            %HR_phase_sample_collection = [HR_phase_sample_collection, phase_sample_store];
            
            % [mean, std] of the estimated parameter
            phase_HR = [((24 / (2*pi)) * phase_mean) * (60 / bin_size), ((24 / (2*pi)) * phase_std) * (60 / bin_size)]; % phase of heart rate circadian rhythm
            
            optfit = [phase_HR(1)];
            optfit_std = [phase_HR(2)/(60 / bin_size)];
            
            %% Compute the uncertainty in the phase parameter
            
            phase_est = mod(optfit(1)/(60 / bin_size), 24);
            %phases = [phases; phase_est];
            phase_est_str = datestr(phase_est/24, 'HH:MM');
            
            %all_fits = [all_fits; optfit] ;
            %all_stds = [all_stds; stdfit];
            
            %% Write day's results to the results.csv file
            
            num_points = length(hr_avg(:,1));
            fprintf(results_file, '%s, %f, %s, %f, %i, %i\n', datestr(dmy_hr(i+1)), optfit(1), phase_est_str, optfit_std(1), num_points, daily_step);
        end
        
    end
    
end

end