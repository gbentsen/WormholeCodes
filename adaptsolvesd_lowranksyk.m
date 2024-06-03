n = 2^15*3*5;                 % number of timepoints
sparsen = 2^7*3*5;            % sparsified number of timepoints (for fermion determinant)
gamma = 5.0;
maxiter = 10000;
targeterror = 1e-6;
kmax = 6;
numbetas = 26;

% betas = 2.^(linspace(-3,10,30));
betas = 2.^(linspace(-3,7,numbetas));
% betas = 2.^(linspace(0,7,numbetas));
% betas = [1 10 50];
% betas = 2.^(-3:10);
% betas = 0.1;

for bb = 1:length(betas)
    
    beta = betas(bb);
    
    if bb == 1
        data(bb) = adaptsolvelowranksyk(n,beta,gamma,maxiter,targeterror,sparsen,kmax);
    else
        oldGw = data(bb-1).Gw / (betas(bb-1)/betas(bb));        % !!! Feed in the old higher-temp solution
        data(bb) = adaptsolvelowranksyk(n,beta,gamma,maxiter,targeterror,sparsen,kmax,oldGw);
    end
    
end

energies = [data.energy]';


energydiffs = real(diff(energies))';
meanbetas = mean([betas(1:end-1);betas(2:end)]);

figure
loglog(betas,energies,'-o');
xlabel('\beta');
ylabel('Energy');

S0 = 1/2*log(2) + cumsum(energydiffs.*meanbetas);

figure
semilogx(1./(betas(2:end)),S0,'-o');
xlabel('Temp T');
ylabel('Entropy Density S/N');

figure
plot(1./(betas(2:end)),S0,'-o');
xlabel('Temp T');
ylabel('Entropy Density S/N');


detdata = reshape([data.det],[kmax,numbetas]);
figure
for kk = 2:kmax
    semilogx(kk./(betas(2:end)),1/2*log(2)-(log(real(detdata(kk,1:end-1)))-log(real(detdata(1,1:end-1))))/(2*(1-kk)),'o-');
    hold on
end
xlabel('Temperature T');
ylabel('MI(A:R)');
title(['Low-Rank SYK, \gamma = ' num2str(gamma)]);
ax = gca; 
ax.FontSize = 14;

    
function data = adaptsolvelowranksyk(n,beta,gamma,maxiter,targeterror,sparsen,kmax,oldGw)
    
    figure
    x = 1.0/16.0;       % It is important that this not begin too high (x = 0.5 is too high)
    preverr = 100.0;

    matsub = zeros(2*n,1);
    Gw = zeros(2*n,1);

    for k=1:2:2*n-1
        if k < n
            matsub(k+1,1) = pi*k/beta;
        else
            matsub(k+1,1) = -pi*(2*n-k)/beta;
        end

        if exist('oldGw','var')     % Feed in the old higher-temp solution
            Gw(k+1,1) = oldGw(k+1,1);
        else
            Gw(k+1,1) = 1./(-1i*matsub(k+1,1));
        end
    end

    data.conds = 0;
    data.divs = 0;
    for ii=1:maxiter

        Gt = fft(Gw) / beta;    % Divide by beta is correct! ([Gw] = +1, [Gt] = 0)

        plot(linspace(0,4*pi,2*n),real(Gt));
        xlabel('\theta');
        ylabel('G(\theta)');
        
% -----------------------------------------------

%         % Regular SYK
%         Sigmat = Gt.^3;
%         Sigmaw = ifft(Sigmat) * beta;
        
        Gsqw = ifft(Gt.^2) * beta;      % Multiply by beta is correct! ([Gsqw] = +1, [Gt] = 0)
        
        if Gsqw(1,1) >= 1.0
            disp('WARNING! Condensation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            Gsqw = Gsqw/(Gsqw(1,1) + 1e-4);
            data.conds = data.conds + 1;
        end
%         assert(Gsqw(1,1) < 1.0);
        
        Fw = 1./(1 - Gsqw);
        Ft = fft(Fw) / beta;
        
        Sigmat = 2*gamma*(Gt.*Ft);
        Sigmaw = ifft(Sigmat) * beta;   % Multiply by beta is correct! ([Sigmaw] = -1, [Sigmat] = -2)
        Sigmaw(abs(Sigmaw) < 1e-12) = 0;
        
% -----------------------------------------------


        assert(~any(Sigmaw(1:2:2*n)));   % ensure that all odd components vanish (ensures that Sigmat is exactly odd in (0,2beta)

        Gwnew = (1-x)*Gw + x * (1./(-1i*matsub-Sigmaw));
        Gwnew(~isfinite(Gwnew))=0;

        err = norm(Gwnew - Gw);

        disp(err);

        if (err < targeterror)
            break
        else
            if err > preverr
                x = x/1.3;
                data.divs = data.divs + 1;
            else
                Gw = Gwnew;
                preverr = err;
            end
        end

    end

%     data.energy = -1/4*sum(Gt.^4)*beta/(2*n);
    data.energy = -gamma/(2*beta)*sum(Gsqw.*Fw);
    data.Gw = Gw;
    
    for ff = 1:kmax
        data.det(ff) = calcdet(Sigmat,ff,n,beta,sparsen);
    end
    
    data.Gbeta2 = Gt(n/2,1);
end