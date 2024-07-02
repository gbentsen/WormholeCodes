alphavals = linspace(0.1,1.0,10);

for ii=1:size(alphavals,2)
    alpha = alphavals(1,ii);
    
    filename = ['SYKFlagpoleGnnABalpha' num2str(alpha) 'New.mat'];
    fprintf(['Loading ' filename]);
    load(filename);
    
    entropyppdata(ii,:) = [data.entropypp];
end