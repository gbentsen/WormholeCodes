n = 2^7*3*5;                 % number of timepoints
maxiter = 1000;
targeterror = 1e-8;
beta = 20;
q = 4;

data = adaptsolvesykmatrix(n,q,beta,maxiter,targeterror);

function data = adaptsolvesykmatrix(n,q,beta,maxiter,targeterror)
    
    figure
    x = 0.5;
    preverr = 100.0;
    
    % Create deltamat
    v = [1 -1];
    v = [v, zeros(1,n-2)];
    deltamat = toeplitz(v,[v(1) fliplr(v(2:end))]);
    deltamat(1,n) = 1;  % antiperiodic boundary conditions

    Gnnp = - transpose(inv(deltamat));

    data.divs = 0;
    for ii=1:maxiter

        plot(linspace(0,2*pi,n),Gnnp(1,:));
        xlabel('\theta');
        ylabel('G(\theta)');
        
% -----------------------------------------------

        Sigmannp = Gnnp.^(q-1);

% -----------------------------------------------

        Gnnpnew = (1-x)*Gnnp + x * (-transpose(inv(deltamat-(beta/n)^2 * Sigmannp)));
        
        diffnnp = Gnnpnew - Gnnp;

        err = abs(trace(diffnnp*conj(diffnnp)))/n;

        disp(err);

        if (err < targeterror)
            break
        else
            if err > preverr
                x = x/2;
                data.divs = data.divs + 1;
            else
                Gnnp = Gnnpnew;
                preverr = err;
            end
        end

    end
    
    data.Gnnp = Gnnp;
end