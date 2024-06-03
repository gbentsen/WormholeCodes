function detval = calcdet(Sigmat,k,n,beta,sparsen)
%CALCDET Summary of this function goes here
%   Detailed explanation goes here
    
    assert(mod(n/sparsen,1)==0);

    Sigmat = Sigmat(1:n/sparsen:end,:);
    n = sparsen;
    
    assert(isequal(size(Sigmat),[2*n 1]));
    
    Sigmatmat = zeros(2*n,n);
    for ii = 0:n-1
        Sigmatmat(:,ii+1) = circshift(Sigmat,ii,1);
    end
    Sigmatmat = Sigmatmat(1:n,:);
    
    partialt = calcpartialtmat(k,n);
    detval = det(partialt - (beta/n)^2*Sigmatmat);

end


function partialt = calcpartialtmat(k,n)

    assert(mod(n,k)==0);
    
    e = ones(n/k,1);
    pt = spdiags([-e e],-1:0,n/k,n/k);
    pt(1,n/k) = 1;
    
    partialt = pt;
    
    if k>1
        for ii = 2:k
            partialt = blkdiag(partialt,pt);
        end
    end
end