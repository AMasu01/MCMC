clear all
clc
close all


Target_distribution = [0.1,0.9]';

dopause = 0;

N_states = numel(Target_distribution);

chain= mcmix(N_states); %Creates a random transition matrix
P = chain.P;

      
G = digraph(1,1:N_states,P(1,:));
for i = 2:N_states
G = addedge(G,i,1:N_states,P(i,:));
end


subplot(221)
h = plot(G,'EdgeLabel',G.Edges.Weight,'LineWidth',5*(G.Edges.Weight)+eps,'EdgeColor',[ 0.5843 0.8157 0.9882]);
title([num2str(N_states), ' State Markov Chain'])



for i = 1:N_states
    for j = 1:N_states
            correction(i,j) = min(1,(Target_distribution(j)*P(j,i))/(Target_distribution(i)*P(i,j)));
    end
end

additional_states = sum(correction~=1);

% Initial (self) edges
G_MH = digraph(num2str(1),num2str(1),P(1,1));
for i = 2:N_states
     G_MH = addedge(G_MH,num2str(i),num2str(i),P(i,i));
end

counter = 0;
for i = 1:N_states
    for j = 1:N_states
        if ~ (i == j) % only moves
           if ~(correction(i,j) == 1)   % do we need an additional state?
               counter  = counter +1; % count additional states
                    G_MH = addedge(G_MH,num2str(i),[num2str(i),'->',num2str(j)]    ,P(i,j)); % Add auxilliary edge
                    G_MH = addedge(G_MH,[num2str(i),'->',num2str(j)] ,num2str(j),correction(i,j)); %
                    G_MH = addedge(G_MH,[num2str(i),'->',num2str(j)] ,num2str(i),1-correction(i,j));
                    % Save the auxilliary nodes for coloring
                    aux_edges{counter} = {{num2str(i),[num2str(i),'->',num2str(j)]},{[num2str(i),'->',num2str(j)] ,num2str(j)},{[num2str(i),'->',num2str(j)] ,num2str(i)}};
           else
                 G_MH = addedge(G_MH,num2str(i),num2str(j),P(i,j));
           end
        end
    end
end


subplot(222)
g = plot(G_MH,'EdgeLabel',G_MH.Edges.Weight,'LineWidth',5*(G_MH.Edges.Weight)+eps,'EdgeColor',[ 0.5843 0.8157 0.9882]);
title(['Corrected ',num2str(N_states), ' State Markov Chain'])

N_draw = 1000; % number of draws

act_node = 1;
act_node_old = 1;
node_sample = zeros(1,N_draw);

act_node_MH = 1;
act_node_old_MH = 1;
node_sample_MH = zeros(1,N_draw);

for i = 1:N_draw
    
subplot(221)
         highlight(h,[act_node_old],'NodeColor','b')
         act_node = randsample(1:N_states,1,true,P(act_node_old,:));
        highlight(h,[act_node],'NodeColor','g')
        highlight(h,act_node,act_node_old,'EdgeColor','red')
        drawnow
        
            if dopause
        pause
        end
        highlight(h,act_node,act_node_old,'EdgeColor',[ 0.5843 0.8157 0.9882])
        act_node_old = act_node;

        node_sample(i) = act_node;

subplot(223)
[a,counts] = histc(node_sample(1:i),1:N_states);
bar(1:N_states,a)


subplot(222)

% first step
         highlight(g,[act_node_old_MH],'NodeColor','b') % Start node
        act_node_MH = randsample(1:N_states,1,true,P(act_node_old_MH,:)); % Sample new node
        
        
        % is it a direct path?
        if ~(correction(act_node_old_MH,act_node_MH) == 1)  && act_node_old_MH ~= act_node_MH%  indirect
            highlight(g,[num2str(act_node_old_MH),'->',num2str(act_node_MH)],'NodeColor','y') % jump to the auxilliary node
            highlight(g,num2str(act_node_old_MH),[num2str(act_node_old_MH),'->',num2str(act_node_MH)],'EdgeColor','red') % mark the path
            drawnow
            if dopause
                pause
            end
            highlight(g,[num2str(act_node_old_MH),'->',num2str(act_node_MH)],'NodeColor','b')
            highlight(g,num2str(act_node_old_MH),[num2str(act_node_old_MH),'->',num2str(act_node_MH)],'EdgeColor',[ 0.5843 0.8157 0.9882])
            
        else   % direct
            highlight(g,num2str(act_node_MH),'NodeColor','y')  % new node
            highlight(g,num2str(act_node_old_MH),num2str(act_node_MH),'EdgeColor','red') % mark the path
            drawnow
            if dopause
                pause
            end
            highlight(g,num2str(act_node_old_MH),num2str(act_node_MH),'EdgeColor',[ 0.5843 0.8157 0.9882])
            highlight(g,num2str(act_node_MH),'NodeColor','b')
        end
        

        

         
% second step         
         if ~(correction(act_node_old_MH,act_node_MH) == 1)    && act_node_old_MH ~= act_node_MH  % if we need a correction 
             inds = randsample([1,0],1,true,[correction(act_node_old_MH,act_node_MH),1-correction(act_node_old_MH,act_node_MH)]); % choose if accept or not
             if inds % if accept                
                 highlight(g,[num2str(act_node_MH)],'NodeColor','y')
                 highlight(g,[num2str(act_node_old_MH)','->',num2str(act_node_MH)],num2str(act_node_MH),'EdgeColor','red')
                 drawnow
                 highlight(g,[num2str(act_node_MH)],'NodeColor','b')
                 highlight(g,[num2str(act_node_old_MH),'->',num2str(act_node_MH)],num2str(act_node_MH),'EdgeColor','g')
             else
                 highlight(g,[num2str(act_node_old_MH)],'NodeColor','y')
                 highlight(g,[num2str(act_node_old_MH),'->',num2str(act_node_MH)],num2str(act_node_old_MH),'EdgeColor','red')
                 drawnow
                   highlight(g,[num2str(act_node_old_MH)],'NodeColor','b')
                  highlight(g,[num2str(act_node_old_MH)','->',num2str(act_node_MH)],num2str(act_node_old_MH),'EdgeColor','g')
                  act_node_MH = act_node_old_MH;
             end
         end
         

    

            if dopause
            pause
            end


            act_node_old_MH = act_node_MH;
            node_sample_MH(i) = act_node_MH;

        
        subplot(224)
        [a,counts] = histc(node_sample_MH(1:i),1:N_states);
        bar(1:N_states,a)
        drawnow

end


