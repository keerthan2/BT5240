clear all
clc

edges = dlmread('webkb-wisc.edges');
u = edges(:,1) + 1; % add 1 to convert from 0 based node indexing to 1 based node indexing
v = edges(:,2) + 1; % add 1 to convert from 0 based node indexing to 1 based node indexing
G = graph(u, v);
plot(G)
%% 2.a
degree_G = degree(G); 
figure
histogram(degree_G,'BinMethod','integers');
title("Degree Distribution")
xlabel("Degree")
%% 2.b
A_sparse = adjacency(G);
cc_G = clustering_coefficients(A_sparse); % MATLAB BGL function

figure
scatter(degree_G, cc_G);
title("Clustering coefficient vs Degree")
xlabel("Degree")
ylabel("Clustering coefficient")

%% 2.c

% Creating 100 random graphs with similar degree distribution
avg_node_degree = mean(degree_G)
avg_cc = mean(cc_G)
% Since the social network is power-law network, we start off with a star
% network and then rewire randomly to produce random network of a similar 
% degree distribution.
[A_s xy] = star_graph(N); 
A_s_sp = sparse(A_s);
R = 100;
avg_deg_rand = zeros(R,N);
N = size(A_s,1);
L_ER = zeros(R,1);
avg_cc_ER = zeros(R,1);
for r=[1:R]
    A_random = WattsStrogatz_network(N,[],1,A_s_sp);
    A_random_sparse = sparse(A_random);
    avg_deg_rand(r,:) = degree(graph(A_random_sparse));
    L_ER(r) = characteristic_path_length(A_random_sparse);
    avg_cc_ER(r) = mean(clustering_coefficients(A_random_sparse)); % MATLAB BGL function
end
avg_rand_deg_dist = mean(avg_deg_rand,1);
figure
histogram(avg_rand_deg_dist,'BinMethod','integers');
title("Mean Degree Distribution of 100 random networks")
xlabel("Degree")

figure
histogram(avg_cc_ER)
hold on
xline(avg_cc,'--','LineWidth',2)
xlabel('Average clustering coefficient')
hold off
legend('Distribution of mean CC in ER','mean CC in social network')
title('Comparing average clustering coefficients of random networks to social network')

L_G = characteristic_path_length(A_sparse)