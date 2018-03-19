% Illustration of the Bootstrap Particle Filter

clear all
close all
clc

T = 100;
THETA_sim = [0.95,1.1,0.2];
simulate_data % simulate data from a simple sv model


THETA = [0.95,1.1,0.2];
DATA = Y_sim;

N.p_x = 1000; % Number of particles

weights_x   = ones(N.p_x,1);

X_old = normrnd(0,1,[N.p_x,1]);
llikeli = 0;


for i = 1:size(DATA,2)
    
    % Transition Equation
    X           = THETA(1)*X_old + THETA(2)*normrnd(0,1,[N.p_x,1]);
    
    % Observation weight
    weights_y   = normpdf(DATA(i),0,exp(X/2)*THETA(3)); % Weights in levels
    
    % Log_likelihood
    icr_likeli = (weights_y'*weights_x)/N.p_x; % Weights in levels
    l_incr = log(icr_likeli);
    llikeli = llikeli + log(icr_likeli);
    
    % Weights
    weights_update     = transform_weights(log(weights_y.*weights_x));
    
    % Resample
    index_resample      = my_rndsamp(weights_update',N.p_x);
    X_old                        = X(index_resample);
    weights_x                 = ones(N.p_x,1);  % We resample every step, thus the weights are equal to unity
    
    % Plots
    if ~exist('fig1')
        fig1 = figure;
        f_size = [14 12];
        fig1.Color = 'w';
        set(fig1,'PaperUnits','inches')
        set(fig1,'PaperSize',f_size)
        fig1.Name = ['IRF'];
        fig1.NumberTitle = 'off';
        fig1.Units = 'inches';
        fig1.Position = [0 0 f_size];
    end
    
    subplot(321)
    scatter(i*ones(size(X)),X,1e3*weights_update+1,'filled','k')
    title('X before resample')
    xlabel('time')
    ylabel('X')
    hold on
    scatter(i,X_sim(i),'filled','r')
    
    subplot(322)
    scatter(i*ones(size(X)),X_old,1e1*1,'filled','k')
    title('X after resample')
    xlabel('time')
    ylabel('X')
    hold on
    hold on
    scatter(i,X_sim(i),'filled','r')
    
    drawnow
    
    subplot(323)
    histogram(X,50)
    title('Histogram of X before resample')
    xlim([-4,4])
    
    subplot(325)
    [~,index_sort] = sort(X);
    plot(X(index_sort),weights_update(index_sort),'-')
    title('Weights of X before resample')
    xlim([-4,4])
    
    subplot(324)
    histogram(X_old,50)
    title('Histogram of X after resample')
    xlim([-4,4])
    
    
    subplot(326)
    plot(i,X_sim(i),'ok')
    hold on
    plot(i,mean(X_old),'og')
    drawnow
    title('X (black) and mean of X* (green)')
    
 %   pause
end
