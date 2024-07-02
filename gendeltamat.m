function deltamat = gendeltamat(len,k)
    v = [1 -1];
    v = [v, zeros(1,len-2)];
    deltamatsingle = toeplitz(v,[v(1) fliplr(v(2:end))]);
    deltamatsingle(1,len) = 1;  % antiperiodic boundary conditions in \beta
    
    deltamatcell = repelem({deltamatsingle}, k);
    deltamat = blkdiag(deltamatcell{:});
end