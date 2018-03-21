

 function [llikeli,llikeli_track] = particle_filter(N_p_x,DATA,THETA)
weights_x = ones(N_p_x,1);

X_old     = normrnd(0,1,[N_p_x,1]);
llikeli   = 0;
  llikeli_track = zeros(1,size(DATA,2));

for i = 1:size(DATA,2)
    
    % Transition Equation
    X           = THETA(1)*X_old + THETA(2)*normrnd(0,1,[N_p_x,1]);
    
    % Observation weight
    weights_y   = normpdf(DATA(i),0,exp(X/2)*THETA(3)); % Weights in levels
    
    % Log_likelihood
    icr_likeli = (weights_y'*weights_x)/N_p_x; % Weights in levels
    l_incr = log(icr_likeli);
    llikeli = llikeli + log(icr_likeli);
    
    % Track the increase in Log-likehood for illustration
    llikeli_track(i) = llikeli;
    
    % Weights
    weights_update     = transform_weights(log(weights_y.*weights_x));
    
    % Resample
    index_resample      = my_rndsamp(weights_update',N_p_x);
    X_old               = X(index_resample);
    weights_x           = ones(N_p_x,1);  % We resample every step, thus the weights are equal to unity
    

end
