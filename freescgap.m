W = 1e5;
T = 1/W;
mu = 5;
iter = 1000;

t = 10;
delta = 5;

iw = 1i*(2*(0:W-1)+1)*pi*T;
Fiw = randn(1,W);
Gsigma = ones(1,W) ./ (iw + mu);

Gexact = -(iw / (2*t^2)) .* (sqrt(conj(iw).*iw + (4*t^2+delta^2)*ones(1,W)) ./ sqrt(conj(iw).*iw + (delta^2)*ones(1,W))-ones(1,W));

for ii = 1:iter
    phi = delta*ones(1,W) + t^2*Fiw;
    
    Gsigmamagn2 = conj(Gsigma) .* Gsigma;
    phimagn2 = conj(phi) .* phi;
    
    GFdenom = 1 + phimagn2 .* Gsigmamagn2;
    
    Giw = Gsigma ./ GFdenom;
    
    Gsigma = ones(1,W) ./ (iw + mu - t^2 * Giw);
    
    Fiw = (phi .* Gsigmamagn2) ./ GFdenom;
    
    mean(mean(Giw))
    
end

figure
zpts = iw(1,1:500:end);
uvals = Giw(1,1:500:end);
plot(1:500:W,imag(uvals),'o')
hold on
plot(imag(Giw))
plot(imag(Gexact))