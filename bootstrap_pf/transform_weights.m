function [ Nweights] = transform_weights( log_weights )
%TRANSFORM_WEIGHTS Summary of this function goes here
%   Detailed explanation goes here

precision                      = 10^(-16);
N                              = length(log_weights);

log_weights                    = log_weights - max(log_weights);
drop                           = log(precision/N);
weights                        = exp(log_weights); % Unnormalized weights
weights(log_weights < drop)    = 0;
weights(isnan(weights))        = 0;
Nweights                       = weights/sum(weights);

end

