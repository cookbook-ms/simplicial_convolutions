% Script to create street network example analysis for London street
% networks
% M. Schaub, June 2018

clc;
clear;
close all;
rng(1223)

data = importdata('LondonEdges.csv');
coordinate = importdata('LondonNodes.csv');
x_coord = coordinate.data(:,2);
y_coord = coordinate.data(:,3);
from = str2double(data.textdata(2:end,1))';
to = str2double(data.textdata(2:end,2))';
num_nodes = max([from, to]);

% num_edges = length(from);
% B1 = sparse(from,1:num_edges,-1,num_nodes,num_edges) ...
%     + sparse(to,1:num_edges,1,num_nodes,num_edges);
% some edges are encoded in both directions -- we need to eliminate them
% for this example, so we construct the adjacency matrix first, and from
% there get the incidence matrix.
A = sparse(from,to,1,num_nodes,num_nodes);
A = A+A';
A(A>1) = 1;
[B1, edgelist] = createIncidenceMatrix(A);
from = edgelist(:,2)';
to = edgelist(:,3)';
num_edges= size(B1,2);

% creat the edge-triangle incidence matrix
[B2, triangle_list] = creat_b2(A,edgelist);

% L1 Hodge Laplacian without faces / Edge Laplacian
L1 = B1'*B1;
% projectors;
Pgrad =  pinv(full(B1))*B1;
Pharm =  eye(size(L1)) -Pgrad;

% Line graph (undirected) adjacency)
AL = abs(L1-2*eye(size(L1)));
LL = diag(sum(AL,2)) - AL;

L0 = B1*B1';


%% Create initial flows and noisy realizations
% random harmonic initial condition
fharm = Pharm*randn(length(L1),1);
fharm = fharm/norm(fharm);

% random low-pass line-graph solution
% first create valid RHS, then solve Linegraph equation
pLL = (eye(size(LL))-ones(size(LL))/num_edges)*randn(length(L1),1);
fLL = 1*(LL\pLL);
fLL = fLL/norm(fLL);

% initial flow a mixture of harmonic and 'linegraph' flows.
flows_noise_free = 10*fharm + fLL;
flows_initial = 1*flows_noise_free + 0.8*randn(num_edges,1);
flows_grad = Pgrad*flows_initial;
flows_harm = Pharm*flows_initial;

p1 = figure;
plt_network(B1,flows_initial,x_coord,y_coord);
axis tight


fnorm = norm(flows_noise_free)
finit_error = norm(flows_noise_free-flows_initial)

plot(flows_noise_free);
hold all;
plot(flows_initial);


%%
% filtering with L1 Laplacian -- gets correct split of flows according to 
% 'Kirchhoffs laws'
best = 9999;
for lambda = 1:1:100
flows_filtered_L1 = (eye(size(L1))+lambda*L1)\flows_initial;
error_L1 = norm(flows_filtered_L1-flows_noise_free);
if error_L1 < best
    best = error_L1;
    lambda_best = lambda;
end
end
best_errorL1 = best
chosen_lambda = lambda_best
flows_filtered_L1 = (eye(size(L1))+chosen_lambda*L1)\flows_initial;

% filtering with AL -- drives all the flows to be equal
figure
plot(flows_noise_free);
hold all;
plot(flows_initial)
best = 9999;
for alpha = 0.01:0.01:10
flows_filtered_LL = (eye(size(L1))+alpha*LL)\flows_initial;
error_LL = norm(flows_filtered_LL-flows_noise_free);
if error_LL < best
    best = error_LL;
    alpha_best = alpha;
end
end

best_errorLL = best 
chosen_alpha = alpha_best
flows_filtered_LL = (eye(size(L1))+chosen_alpha*LL)\flows_initial;


best = 9999;
for alpha = 0.01:0.01:10
    for lambda = 1:1:100
        
        flows_filtered_mixed = (eye(size(L1))+alpha*LL+lambda*L1)\flows_initial;
        error_mixed = norm(flows_filtered_mixed-flows_noise_free);
        if error_mixed < best
            best = error_mixed;
            alpha_best = alpha;
            lambda_best = lambda;
        end
    end
end
best
alpha_best
lambda_best
flows_filtered_mixed = (eye(size(L1))+alpha_best*LL+lambda_best*L1)\flows_initial;

%% smoothing 
figure;
plot(flows_noise_free);
hold all;
plot(flows_initial);
hold on;
flows = (eye(size(L1))-0.2*L1)*flows_initial;
for t = 1:100
    flows = (eye(size(L1))-0.2*L1)*flows;
    plot(flows); hold on; 
    e(t) = flows'*L1*flows;
end

%% Write out results to file
M = [from', to', flows_noise_free, flows_initial, flows_filtered_L1, flows_filtered_LL,flows_filtered_mixed,flows_grad,flows_harm];
cHeader = {'#NoiseFreeFlow','NoisyFlow','FilteredFlowL1','FilteredFlowLL','Filtered_mixed','FlowGrad','FlowHarm'};
textHeader = strjoin(cHeader, ',');
fid = fopen('StreetNetworkLondon.csv','w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('StreetNetworkLondon.csv',M,'-append','delimiter',',')
