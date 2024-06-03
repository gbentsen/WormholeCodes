function actionbyn = calcactionbyn(GnnpA, GnnpB, alpha, beta, n, k)

    betabynsqr = (beta/n)^2;

    % Create deltamatA
    v = [1 -1];
    v = [v, zeros(1,k*n-2)];
    deltamatA = toeplitz(v,[v(1) fliplr(v(2:end))]);
    deltamatA(1,k*n) = 1;  % antiperiodic boundary conditions in 2\beta

    % Create deltamatB
    v = [1 -1];
    v = [v, zeros(1,n-2)];
    deltamatB = toeplitz(v,[v(1) fliplr(v(2:end))]);
    deltamatB(1,n) = 1;  % antiperiodic boundary conditions in \beta
    deltamatB = blkdiag(deltamatB,deltamatB);


    Sigmannp = ( alpha * GnnpA + (1-alpha) * GnnpB ).^3;

    Knnp = ( alpha * GnnpA + (1-alpha) * GnnpB ).^2;

    detA = det(deltamatA - betabynsqr * Sigmannp);
    detB = det(deltamatB - betabynsqr * Sigmannp);
    energytrace = betabynsqr * trace(alpha * transpose(Sigmannp)*GnnpA + (1-alpha) * transpose(Sigmannp)*GnnpB - 1/4 * transpose(Knnp)*Knnp);


    actionbyn = -1/2*alpha * log(detA) - 1/2*(1-alpha) * log(detB) + 1/2*energytrace;
        
end