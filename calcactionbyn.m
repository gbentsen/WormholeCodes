function actionbyn = calcactionbyn(GnnpA, GnnpB, alpha, beta, n, k, kA, kB)

    if (~exist('kA','var'))
        kA = 1;
    end
    if (~exist('kB','var'))
        kB = 2;
    end
    if (mod(k,kA) ~= 0 || mod(k,kB) ~= 0)
        error('k must be divisible by kA and kB');
    end
    
    betabynsqr = (beta/n)^2;

    % Create deltamatA
    deltamatA = gendeltamat(k/kA*n,kA);
    
    % Create deltamatB
    deltamatB = gendeltamat(k/kB*n,kB);


    Sigmannp = ( alpha * GnnpA + (1-alpha) * GnnpB ).^3;

    Knnp = ( alpha * GnnpA + (1-alpha) * GnnpB ).^2;

    detA = det(deltamatA - betabynsqr * Sigmannp);
    detB = det(deltamatB - betabynsqr * Sigmannp);
    energytrace = betabynsqr * trace(alpha * transpose(Sigmannp)*GnnpA + (1-alpha) * transpose(Sigmannp)*GnnpB - 1/4 * transpose(Knnp)*Knnp);


    actionbyn = -1/2*alpha * log(detA) - 1/2*(1-alpha) * log(detB) + 1/2*energytrace;
        
end