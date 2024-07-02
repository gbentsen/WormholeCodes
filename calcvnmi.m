kmax = 6;
numbetas = 30;

% detdata = reshape([data.det],[kmax,numbetas]);

% m2 = (1/2*log(2)-(log(real(detdata(2,:)))-log(real(detdata(1,:))))/(2*(1-2)));

usqr = exp(m2)-1;
mVNpred = log(2) + 1/2*(1+sqrt(usqr)).*log(1/2*(1+sqrt(usqr))) + 1/2*(1-sqrt(usqr)).*log(1/2*(1-sqrt(usqr)));