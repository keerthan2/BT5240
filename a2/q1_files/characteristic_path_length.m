function L = characteristic_path_length(A)
    % A should be sparse
    N = size(A,1);
    d_array = zeros(N,1);
    for node=[1:N]
        d_ij = shortest_paths(A,node); % MATLAB BGL function
        d_ij = d_ij(d_ij ~= 0);
        d_ij = d_ij(~isinf(d_ij));
        if numel(d_ij)
            d_array(node) = mean(d_ij);
        end
    end
    L = mean(d_array);
end