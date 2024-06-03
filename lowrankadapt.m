maxiter = 10000;
targeterror = 1e-6;
x = 1e-2;
gamma = 1;

N = 2^22;
% betas = [0.125 0.25 0.5 1 2 4 8 16 32 64 128 256];
betas = 2.^(linspace(-3,6,20));

% betas = 50;

energies = zeros(length(betas),1);

for bb = 1:length(betas)
    beta = betas(bb);
    
    phase = exp(1i*pi*(0:(N-1))/N);
    iw = 1i*(2*pi/beta)*((-floor(N/2):N-1-floor(N/2))+1/2);
%     iw = 1i*(2*pi/beta)*((-floor(N/2)+1:N-1-floor(N/2)+1)+1/2);
%     iw = 1i*(2*pi/beta)*((-ceil(N/2)+1:N-ceil(N/2))+1/2);
%     iw = 1i*(2*pi/beta)*((0:N-1)+1/2);

    Gw = 1./(-iw);

    preverr = 10;

    for ii = 1:maxiter

        Gt = fft(Gw) .* conj(phase) / beta;
        
        Gsqw = ifft(Gt.^2) * beta;
        Fw = 1./(1 - Gsqw);
        Ft = fft(Fw) / beta;
        
        Sigmat = 2*gamma*(Gt.*Ft);
        Sigmaw = ifft(Sigmat .* phase) * beta;
        Gwnew = (1-x)*Gw + x * (1./(-iw-Sigmaw));

        err = norm(Gwnew - Gw);

        disp(err);

        if (err < targeterror)
            break
        else
            if err > preverr
                x = x/2;
            else
                Gw = Gwnew;
                preverr = err;
            end
        end
    end

    energies(bb) = -1/2*gamma*sum(Gsqw.*Fw);
end

energydiffs = real(diff(energies))';
meanbetas = mean([betas(1:end-1);betas(2:end)]);

figure
loglog(betas,energies,'-o');

S0 = 1/2*log(2) + cumsum(energydiffs.*meanbetas);