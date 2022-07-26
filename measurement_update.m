function [x_corrected, M_corrected] = measurement_update(measurement_inst, xM_phase)
        
    measurement = measurement_inst.hr_phase;
    measurement_fun = @measurement_function;
    measurement_noise = measurement_inst.k_observe;
    
%   Input Arguments:
%   xS0: average and a square root of predicted covariance matrix of class
%   mean_covariance_sqrt_cls
%   measurement: measurement at step
%   measurement_fun: (non-stochastic) measurement function z = f(x)
%   measurement_noise: Covariance matrix of measurement noise (not a square
%   root)
%   Output Argument:
%   xS_corrected: corrected average and square root covariance matrix of
%   class mean_covariance_sqrt_cls
    
   [x_corrected, M_corrected] = square_root_cubature_measurement_update(xM_phase,measurement,measurement_fun,measurement_noise);
 
end