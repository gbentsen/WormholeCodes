n = 2^17*3*5;                 % number of timepoints
sparsen = 2^7*3*5;            % sparsified number of timepoints (for fermion determinant)
maxiter = 1000;
targeterror = 1e-8;
kmax = 6;
numbetas = 1;
q = 4;

% betas = 5.2029*2;
betas = 2^(-3);

% betas = 2.^(linspace(-3,10,numbetas));
% betas = 2.^(-3:10);
% betas = 0.1;
% betas = 0.1;

for bb = 1:length(betas)
    
    beta = betas(bb);
    
    data(bb) = adaptsolvesyk(n,q,beta,maxiter,targeterror,sparsen,kmax);
    
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
title('SYK');
ax = gca; 
ax.FontSize = 14;


m2 = 1/2*log(2)-(log(real(detdata(2,:)))-log(real(detdata(1,:))))/(2*(1-2));
usqr = exp(m2)-1;
mVNpred = log(2) + 1/2*(1+sqrt(usqr)).*log(1/2*(1+sqrt(usqr))) + 1/2*(1-sqrt(usqr)).*log(1/2*(1-sqrt(usqr)));
figure
loglog(2./(betas(13:end)),mVNpred(13:end),'o-')
figure
plot(log(2./(betas(15:end))),log(mVNpred(15:end)),'o-')


function data = adaptsolvesyk(n,q,beta,maxiter,targeterror,sparsen,kmax)
    
    figure
    x = 0.5;
    preverr = 100.0;

    matsub = zeros(2*n,1);
    Gw = zeros(2*n,1);

    for k=1:2:2*n-1
        if k < n
            matsub(k+1,1) = pi*k/beta;
        else
            matsub(k+1,1) = -pi*(2*n-k)/beta;
        end

        Gw(k+1,1) = 1./(-1i*matsub(k+1,1));
    end

    data.divs = 0;
    for ii=1:maxiter

        Gt = fft(Gw) / beta;    % Divide by beta is correct! ([Gw] = +1, [Gt] = 0)

        plot(linspace(0,4*pi,2*n),real(Gt));
        xlabel('\theta');
        ylabel('G(\theta)');
        
% -----------------------------------------------

        Sigmat = Gt.^(q-1);
        Sigmaw = ifft(Sigmat) * beta;   % Multiply by beta is correct! ([Sigmaw] = -1, [Sigmat] = -2)

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
                x = x/2;
                data.divs = data.divs + 1;
            else
                Gw = Gwnew;
                preverr = err;
            end
        end

    end

    data.energy = -1/q*sum(Gt.^q)*beta/(2*n);   % dt = beta/n, divide by 2 because Gt actually goes from (0,2beta)
    
    for ff = 1:kmax
        data.det(ff) = calcdet(Sigmat,ff,n,beta,sparsen);
    end
    
    data.Gbeta2 = Gt(n/2,1);
end