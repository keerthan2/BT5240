function A = WattsStrogatz_network(N,k,p,A)
    if nargin < 4
        A = regular_lattice(N,k);
    end
    for node=[1:N]
        edges = find(A(node,:) == 1);
        non_edges = find(A(node,:) == 0);
        non_edges = non_edges(non_edges ~= node);
        n_edges = numel(edges);
        n_non_edges = numel(non_edges);
        isRewire = rand(n_edges,1) < p;
        edges_to_rewire = edges(edges.*isRewire' ~= 0);
        n_etr = numel(edges_to_rewire);
        idx = randperm(n_non_edges,n_etr);
        new_targets = non_edges(idx);
        A(node,new_targets) = 1;
        A(new_targets,node) = 1;
        A(node,edges_to_rewire) = 0;
        A(edges_to_rewire,node) = 0; 
    end
end

function A = regular_lattice(N,k)
    A = zeros(N);
    for node=1:N
        high = [];
        low = [];
        for j=1:k
            high = [high mod(node+j,N)];
            low = [low mod(node-j,N)];
        end
        high(high==0) = N;
        low(low==0) = N;
        A(node,high) = 1;
        A(node,low) = 1;
    end
end