function data = adaptsolvesykmatrixflagpole(n,q,betas,maxiter,x,targeterror,alphas,k,calcaction,kA,kB)
    
    for aa = 1:length(alphas)
        for bb = 1:length(betas)

            alpha = alphas(aa);
            beta = betas(bb);

            fprintf('alpha = %f, beta = %f\n',alpha,beta);
    %         if (bb == 1)
                data(aa,bb) = adaptsolvesykmatrixflagpolesingle(n,q,beta,maxiter,x,targeterror,alpha,k,calcaction,kA,kB);
    %         else
    %             oldGnnpA = data(aa,bb-1).GnnpA;        % !!! Feed in the old higher-temp solution
    %             oldGnnpB = data(aa,bb-1).GnnpB;
    %             data(aa,bb) = adaptsolvesykmatrixflagpole(n,q,beta,maxiter,targeterror,alpha,k,oldGnnpA,oldGnnpB);
    %         end

        end 
    end
end

function data = adaptsolvesykmatrixflagpolesingle(n,q,beta,maxiter,x,targeterror,alpha,k,calcaction,kA,kB,oldGnnpA,oldGnnpB)
    
    if (~exist('kA','var'))
        kA = 1;
    end
    if (~exist('kB','var'))
        kB = 2;
    end
    if (mod(k,kA) ~= 0 || mod(k,kB) ~= 0)
        error('k must be divisible by kA and kB');
    end
    
    preverr = 1000000.0;
    
    betabynsqr = (beta/n)^2;
    
    % Create deltamatA
    deltamatA = gendeltamat(k/kA*n,kA);
    
    % Create deltamatB
    deltamatB = gendeltamat(k/kB*n,kB);

    
    if (exist('oldGnnpA','var') && exist('oldGnnpB','var'))     % Feed in the old higher-temp solution
        GnnpA = oldGnnpA;
        GnnpB = oldGnnpB;
    else
%         GnnpA = - transpose(inv(deltamatA)) + (1/2*rand(k*n,k*n)-1/4);
%         GnnpB = - transpose(inv(deltamatB)) + (1/2*rand(k*n,k*n)-1/4);
        GnnpA = - transpose(inv(deltamatA));
        GnnpB = - transpose(inv(deltamatB));
%         GnnpA = rand(k*n,k*n)-1/2;
%         GnnpB = rand(k*n,k*n)-1/2;
    end

    data.divs = 0;
    for ii=1:maxiter

%         if preverr < 0.1
%             x = 0.25/8;
%         end
            
%         plot(linspace(0,2*pi,k*n),GnnpA(1,:),'o-');
%         xlabel('\theta');
%         ylabel('G(\theta)');
        
% -----------------------------------------------

        Sigmannp = ( alpha * GnnpA + (1-alpha) * GnnpB ).^(q-1);

% -----------------------------------------------

        propA = deltamatA - betabynsqr * Sigmannp;
        propB = deltamatB - betabynsqr * Sigmannp;
        
        ntrinvpropA = -transpose(inv(propA));
        ntrinvpropB = -transpose(inv(propB));

        GnnpAnew = (1-x)*GnnpA + x * (ntrinvpropA);
        GnnpBnew = (1-x)*GnnpB + x * (ntrinvpropB);
        
%         GnnpBnew(1:n,n+1:2*n) = 0;
%         GnnpBnew(n+1:2*n,1:n) = 0;
        
        diffnnpA = GnnpAnew - GnnpA;
        diffnnpB = GnnpBnew - GnnpB;

        errA = abs(trace(diffnnpA*diffnnpA'))/n;
        errB = abs(trace(diffnnpB*diffnnpB'))/n;
        
        err = max(errA,errB);
        
        
        absdiffnnpA = GnnpAnew - ntrinvpropA;
        absdiffnnpB = GnnpBnew - ntrinvpropB;

        abserrA = abs(trace(absdiffnnpA*absdiffnnpA'))/n;
        abserrB = abs(trace(absdiffnnpB*absdiffnnpB'))/n;
        
        abserr = max(abserrA,abserrB);
        
        data.abserr = abserr;
        
        if calcaction
            actionbyn = calcactionbyn(GnnpAnew, GnnpBnew, alpha, beta, n, k);

            disp([x err abserr actionbyn]);
        else
            disp([x err abserr]);
        end

        if (err < targeterror)
            break
        else
            if err > preverr
%                 x = x/2;
                data.divs = data.divs + 1;
                break
            else
                GnnpA = GnnpAnew;
                GnnpB = GnnpBnew;
                preverr = err;
            end
        end

    end
    
    data.GnnpA = GnnpA;
    data.GnnpB = GnnpB;
end