% Standard Stochastic Volatility Model

X_sim = 0;
Y_sim = exp(X_sim/2)*THETA_sim(3)*normrnd(0,1);


for t = 2:T
X_sim(t)           = THETA_sim(1)*X_sim(t-1) + THETA_sim(2)*normrnd(0,1,[1,1]);
Y_sim(t)           = exp(X_sim(t)/2)*THETA_sim(3)*normrnd(0,1,[1,1]);

end

% figure
% plot(1:T,X_sim,'-k',1:T,Y_sim,'-.g')
% legend('X','Y')