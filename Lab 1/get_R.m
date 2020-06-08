function [ R ] = get_R(mic_var)
% Returns the covariance matrix R
    nr_of_mics = length(mic_var);
    
    % Let all elements in matrix R be the variance of microphone 8
    R = mic_var(7)*ones(6,6);

    i = 1;
    for k = 1:nr_of_mics-1
       % Update the values for the diagonal in matrix R
       y_var =  mic_var(7) + mic_var(k);
       R(i,k) = y_var;
       i = i+1;
    end
    
end