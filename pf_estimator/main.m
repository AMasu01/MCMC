% Illustration of the Bootstrap Particle Filter

clear all
close all
clc

% Simulate data
T = 100;
THETA_sim = [0.95,1.1,0.2];
simulate_data % simulate data from a simple sv model
DATA      = Y_sim;


%% 1. Illustrate the variance of the likelihood estimator
THETA     = [0.95,1.1,0.2];

N_p_x     = [100,200,500]; % Number of particles
N_variants =  numel(N_p_x);
N_run       = 100; % Number of runs
loglikelihood_tracked = zeros(N_variants,N_run);

for j = 1:N_variants
    for i = 1:N_run
    loglikelihood_tracked(j,i) = particle_filter(N_p_x(j),DATA,THETA);
    end
end
  if ~exist('fig1')
        fig1 = figure;
        f_size = [14 12];
        fig1.Color = 'w';
        set(fig1,'PaperUnits','inches')
        set(fig1,'PaperSize',f_size)
        fig1.Name = ['PF_illustration'];
        fig1.NumberTitle = 'off';
        fig1.Units = 'inches';
        fig1.Position = [0 0 f_size];
  end
  
  for j = 1:N_variants
      [a,b] = ksdensity(loglikelihood_tracked(j,:));
    plot(b,a,'DisplayName',[num2str(N_p_x(j)),' Particles',])
    legend('-DynamicLegend')
    hold on
  end
  title('Log-Likelihood Estimator')
  hold off

 
  pause
  
  
  %% 2. Check cuts through the loglikelihood surface
  N_run                              = 50;
  N_theta                           = numel(THETA);
  names                             = {'\rho','\sigma_y','\sigma_x'};
  loglikelihood_tracked    = zeros(N_variants,N_run,size(DATA,2));
  
  for j = 1:N_theta
      theta_vec(j,:) = linspace(0.5*THETA(j),1.5*THETA(j),N_run);
      for i = 1:N_run
          Theta_use                                      =     THETA;
          Theta_use(j)                                   = theta_vec(j,i);
          [~,loglikelihood_tracked(j,i,:)]   = particle_filter(2000,DATA,Theta_use);
      end
      
      
  end
  
  % Plots
  for k = 1:size(DATA,2)
      for j = 1:N_theta
          subplot(3,1,j)
          plot(theta_vec(j,:) ,squeeze(loglikelihood_tracked(j,:,k)),'-b')
          title(['Likelihood Cut: ',names{j},' at t = ', num2str(k)])
      end
      clc
      fprintf('PRESS ANY KEY')
      pause
      clc
  end

    