% Illustrates Importance Sampling using a bivariate distribution
clc 
clear all
close all

 if ~exist('fig1')
        fig1 = figure;
        f_size = [14 12];
        fig1.Color = 'w';
        set(fig1,'PaperUnits','inches')
        set(fig1,'PaperSize',f_size)
        fig1.Name = ['Importance Sampling'];
        fig1.NumberTitle = 'off';
        fig1.Units = 'inches';
        fig1.Position = [0 0 f_size];
 end
  
 
 
target_distribution = @(x) [gampdf(x(:,1),1,1),normpdf(x(:,2),0,0.2)];
importance_distribution = @(x) [normpdf(x(:,1),0,1),normpdf(x(:,2),0,1)];


% Draw from the importance distribution
N = 1000; % Choose the numer of sampled points

draw_id = [normrnd(0,1,[N,1]),normrnd(0,1,[N,1])];

subplot(321)
    scatter(draw_id(:,1),draw_id(:,2),10e2*1/N,'r')
    title('Importance Distribution & Sample Points')

     x_ax = linspace(min(draw_id(:,1))*0.9,max(draw_id(:,1))*1.1,100)  ;  %// x axis
    y_ax =  linspace(min(draw_id(:,2))*0.9,max(draw_id(:,2))*1.1,100); %// y axis   
    [X ,Y] = meshgrid(x_ax,y_ax); %// all combinations of x, y
    XY = [X(:),Y(:)];
    Z = prod(importance_distribution(XY),2);
    Z = reshape(Z,size(X));
    hold on
    contour(X,Y,Z,20)

subplot(322)
    hist3(draw_id,[20,20])
    [~,edg] = hist3(draw_id,[20,20]);
    title('Histogram Importance Sample')

subplot(323)
    Z = prod(target_distribution(XY),2);
    Z = reshape(Z,size(X));
    contour(X,Y,Z,20)
    hold on
    importance_weights = prod(target_distribution(draw_id),2)./prod(importance_distribution(draw_id),2);
    scatter(draw_id(:,1),draw_id(:,2),10e2*importance_weights/sum(importance_weights)+eps,'r')
    title('Target Distribution & Weighted Sample Points')

subplot(324)
    mesh(X,Y,Z)
    title('Target Distribution')
    xlim([min(edg{1}),max(edg{1})])
    ylim([min(edg{2}),max(edg{2})])

subplot(325)
    contour(X,Y,Z,20)
    hold on
    title('Target Distribution & Resampled Points')
    indices_resample = randsample(1:N,N,true,importance_weights/sum(importance_weights));
    scatter(draw_id(indices_resample,1),draw_id(indices_resample,2),10e2*1/N,'r')

subplot(326)
    hist3(draw_id(indices_resample,:),[20,20]);
    title('Histogramm Target Sample')
    xlim([min(edg{1}),max(edg{1})])
    ylim([min(edg{2}),max(edg{2})])
