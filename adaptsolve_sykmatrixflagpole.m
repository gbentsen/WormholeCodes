n = 2^7*3*5;     % number of timepoints
% n = 2^11;
maxiter = 100;
targeterror = 1e-5;
% beta = 1;
q = 4;
k = 2;

kA = 1;
kB = 2;

% alpha = 0.1;

x = 0.25;

calcaction = false;

% numalphas = 11;
% alphas = linspace(0,1,numalphas);
alphas = [0.99];

numbetas = 30;
betas = 2.^(linspace(-3,10,numbetas));
% betas = [betas(22)];
% numbetas = 1;
% betas = [0.1 1 10 100];
betas = betas(10);

% alphas = [0.5];
% betas = [1];

data = adaptsolvesykmatrixflagpole(n,q,betas,maxiter,x,targeterror,alphas,k,calcaction,kA,kB);