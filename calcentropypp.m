for aa = 1:length(alphas)
    for bb = 1:length(betas)
        
        GnnpAalpha0 = dataalpha0(aa,bb).GnnpA;
        GnnpBalpha0 = dataalpha0(aa,bb).GnnpB;
        
        GnnpA = data(aa,bb).GnnpA;
        GnnpB = data(aa,bb).GnnpB;
        
        alpha = alphas(aa);
        beta = betas(bb);
    
        actionbynalpha0 = calcactionbyn(GnnpAalpha0, GnnpBalpha0, 0.0, beta, n, k);
        actionbyn = calcactionbyn(GnnpA, GnnpB, alpha, beta, n, k);
        
        data(aa,bb).entropypp = actionbyn - actionbynalpha0;
        
        disp(data(aa,bb).entropypp);
        
    end 
end

