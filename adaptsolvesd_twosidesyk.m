n = 2^20;                   % number of timepoints
sparsen = 2^10;              % timepoints to save in data
mu = 0.5;
maxiter = 10000;
targeterror = 1e-8;
% numbetas = 30;
nummus = 30;
% numxs = 30;

% couplingx = 0.011;
% J2 = sqrt(couplingx);
% J4 = sqrt(1-couplingx);
% couplingx = 0.0;
J2 = 0.0;
J4 = 1.0;
% couplings = linspace(0,1,numxs);

beta = 2^10;
% betas = 2.^(linspace(-3,10,numbetas));
% betas = 2.^(-3:10);
% betas = 0.1;
% betas = 50;

mus = linspace(0,0.01,nummus);

tic

for bb = 1:length(mus)
    
%     beta = betas(bb);
    mu = mus(bb);

%     couplingx = couplings(bb);
%     J2 = sqrt(couplingx);
%     J4 = sqrt(1-couplingx);
    
    data(bb) = adaptsolvesyk(n,beta,mu,J2,J4,maxiter,targeterror,sparsen);
    
end

toc

    
function data = adaptsolvesyk(n,beta,mu,J2,J4,maxiter,targeterror,sparsen)
    
    figure
    x = 0.3;
    preverr = 100.0;

    matsub = zeros(2*n,1);
    GLLw = zeros(2*n,1);
    GRRw = zeros(2*n,1);
    GLRw = zeros(2*n,1);
    GRLw = zeros(2*n,1);

    for k=1:2:2*n-1
        if k < n
            matsub(k+1,1) = pi*k/beta;
        else
            matsub(k+1,1) = -pi*(2*n-k)/beta;
        end

        GLLw(k+1,1) = 1./(-1i*matsub(k+1,1));
        GRRw(k+1,1) = 1./(-1i*matsub(k+1,1));
        GLRw(k+1,1) = 0;
        GRLw(k+1,1) = 0;
    end

    data.divs = 0;
    for ii=1:maxiter

        GLLt = fft(GLLw) / beta;    % Divide by beta is correct! ([Gw] = +1, [Gt] = 0)
        GRRt = fft(GRRw) / beta;
        GLRt = fft(GLRw) / beta;
        GRLt = fft(GRLw) / beta;

        plot(linspace(0,4*pi,2*n),real(GLLt));
        xlabel('\theta');
        ylabel('GLL(\theta)');
        
        plot(linspace(0,4*pi,2*n),imag(GLRt));
        xlabel('\theta');
        ylabel('GLR(\theta)');
        
% -----------------------------------------------

        SigmaLLt = J4*2*GLLt.^3 + J2*2*GLLt;
        SigmaRRt = J4*2*GRRt.^3 + J2*2*GRRt;
        SigmaLRt = J4*2*GLRt.^3 - J2*2*GLRt;    % Factor of (-1)^{q/2} because two sides evolve oppositely
        SigmaRLt = J4*2*GRLt.^3 - J2*2*GRLt;
        
%         % Pure J4
%         SigmaLLt = J4*2*GLLt.^3;
%         SigmaRRt = J4*2*GRRt.^3;
%         SigmaLRt = J4*2*GLRt.^3;
%         SigmaRLt = J4*2*GRLt.^3;
        
        SigmaLLw = ifft(SigmaLLt) * beta;   % Multiply by beta is correct! ([Sigmaw] = -1, [Sigmat] = -2)
        SigmaRRw = ifft(SigmaRRt) * beta;
        SigmaLRw = ifft(SigmaLRt) * beta - (1i * mu * repmat([0; 1],n,1));   % Maldacena-Qi coupling
        SigmaRLw = ifft(SigmaRLt) * beta + (1i * mu * repmat([0; 1],n,1));
        
% -----------------------------------------------


        assert(~any(SigmaLLw(1:2:2*n)));   % ensure that all odd components vanish (ensures that Sigmat is exactly odd in (0,2beta)
        assert(~any(SigmaRRw(1:2:2*n)));
        assert(~any(SigmaLRw(1:2:2*n)));
        assert(~any(SigmaRLw(1:2:2*n)));
        
        denom = (-1i*matsub-SigmaLLw).*(-1i*matsub-SigmaRRw) - SigmaLRw .* SigmaRLw;
        
        GLLwnew = (1-x)*GLLw + x * (-1i*matsub-SigmaRRw) ./ denom;
        GRRwnew = (1-x)*GRRw + x * (-1i*matsub-SigmaLLw) ./ denom;
        GLRwnew = (1-x)*GLRw + x * SigmaLRw ./ denom;
        GRLwnew = (1-x)*GRLw + x * SigmaRLw ./ denom;
        
        GLLwnew(~isfinite(GLLwnew))=0;
        GRRwnew(~isfinite(GRRwnew))=0;
        GLRwnew(~isfinite(GLRwnew))=0;
        GRLwnew(~isfinite(GRLwnew))=0;

        err = norm(GLLwnew - GLLw) + norm(GRRwnew - GRRw) + norm(GLRwnew - GLRw) + norm(GRLwnew - GRLw);

        disp(err);

        if (err < targeterror)
            break
        else
            if err > preverr
                x = x/1.3;
                data.divs = data.divs + 1;
            else
                GLLw = GLLwnew;
                GRRw = GRRwnew;
                GLRw = GLRwnew;
                GRLw = GRLwnew;
                preverr = err;
            end
        end

    end
    
    assert(mod(n/sparsen,1)==0);
    data.GLLt = GLLt(1:n/sparsen:end,:);
    data.GRRt = GRRt(1:n/sparsen:end,:);
    data.GLRt = GLRt(1:n/sparsen:end,:);
    data.GRLt = GRLt(1:n/sparsen:end,:);
    
    data.energyLL = -1/4*sum(GLLt.^4)*beta/(2*n) - 1/4*sum(GLRt.^4)*beta/(2*n);   % dt = beta/n, divide by 2 because Gt actually goes from (0,2beta)
    
    
%     data.energy = 
    
end