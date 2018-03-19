% This code illustrates the problem of "uninformative" priors and model
% implied conditions.

% Bivariate prior
% prior_x = @(x) unifrnd(-1,1,[x,1]);
% prior_y = @(x) unifrnd(0,1,[x,1]);

%Alternative prior, uncomment to use
prior_x = @(x) normrnd(0,1,[x,1]);
prior_y = @(x) normrnd(0,1,[x,1]);

% Sample
N_sample = 50000;

x_sample_prior = prior_x(N_sample);
y_sample_prior = prior_y(N_sample);
joint_sample_prior =  [x_sample_prior,y_sample_prior];

condition = y_sample_prior > abs(x_sample_prior.^0.3);  % Model dependent

joint_sample_prior_condition = joint_sample_prior(condition,:);

fig1 = figure
f_size = [12 7];
fig1.Color = 'w';
set(fig1,'PaperUnits','inches')
set(fig1,'PaperSize',f_size)
fig1.Name = ['IRF'];
fig1.NumberTitle = 'off';
fig1.Units = 'inches';
fig1.Position = [0 0 f_size];
    
    
subplot(321)
histogram(x_sample_prior,50)
title('p(X)')
subplot(323)
histogram(y_sample_prior,50)
title('p(Y)')
subplot(325)
hist3(joint_sample_prior,[50,50])
title('p(X,Y)')

subplot(322)
histogram(joint_sample_prior_condition(:,1),50)
title('p(X) with condition')
subplot(324)
histogram(joint_sample_prior_condition(:,2),50)
title('p(Y) with condition')
subplot(326)
hist3(joint_sample_prior_condition,[50,50])
title('p(X,Y) with condition')
