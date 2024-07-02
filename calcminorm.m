for ii=1:9
    minorm(ii,:) = (entropyppdata(ii,:) + entropyppdata(10,:) - entropyppdata(10-ii,:))/(ii/10.0);
end

for ii=1:9
    miunnorm(ii,:) = (entropyppdata(ii,:) + entropyppdata(10,:) - entropyppdata(10-ii,:));
end