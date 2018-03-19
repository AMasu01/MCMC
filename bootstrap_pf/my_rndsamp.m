function [ y] = my_rndsamp(weight,k,U)
%MY_RNDSAMP samples indices based on weight, this the lean version of
%MATLABS "randsample"

if nargin == 2
U = rand(k,1);
end

edges = min([0 cumsum(weight)],1); % protect against accumulated round-off
edges(end) = 1; % get the upper edge exact
%[~, y] = histc(U,edges); % SLOW!
y = discretize(U, edges); % FAST!
end

