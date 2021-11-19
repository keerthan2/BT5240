clc
clear all

N = 100;
R = 100; % Number of realizations of ER network
k = 5;
p = 0.7;
rng(0)
A = WattsStrogatz_network(N,k,0.7);
A_sparse = sparse(A);
L_WS = characteristic_path_length(A_sparse)
avg_cc_WS = mean(clustering_coefficients(A_sparse)) % MATLAB BGL function
G2 = graph(A);
plot(G2);

% generate 100 random networks and compute parameters
L_ER = zeros(R,1);
avg_cc_ER = zeros(R,1);
for r=[1:R]
    A_random = WattsStrogatz_network(N,k,1);
    A_random_sparse = sparse(A_random);
    L_ER(r) = characteristic_path_length(A_random_sparse);
    avg_cc_ER(r) = mean(clustering_coefficients(A_random_sparse)); % MATLAB BGL function
end

figure
subplot(1,2,1)
histogram(avg_cc_ER)
hold on
xline(avg_cc_WS,'--','LineWidth',2)
xlabel('Average clustering coefficient')
hold off
legend('Distribution of mean CC in ER','mean CC in rewired network')
title('Comparing average clustering coefficients of random networks to rewired network')


subplot(1,2,2)
histogram(L_ER)
hold on
xline(L_WS,'--','LineWidth',2)
xlabel('Characteristic path length')
hold off
legend('Distribution of CL in ER','CL in rewired network')
title('Comparing characteristic path length of random networks to rewired network')
