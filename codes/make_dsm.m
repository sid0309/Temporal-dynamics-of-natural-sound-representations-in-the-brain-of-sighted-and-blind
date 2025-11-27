function [dsm, avg_dsm,comparison] = make_dsm(data,n_cat,factor)
% Written by Siddharth Talwar
% Last edited on 26-11-2025
% makes dissimilarity matrices. Input data should be 2dimensional -
% comparison x time. Also computes avg of dsm at each time point using a 
% moving window. The length of moving window is given by factor which is 
% number of neighbouring/ adjacent bins on 1 side of the timebin.

if length(size(data)) > 2
    error('the input should be 2 dim: comparisons x time');
end

comparison = flip(combnk(1:n_cat,2));
if comparison(1,1)~=1 % to make sure the combinations start with the 1st
    comparison = flip(comparison);
end

dsm = zeros(n_cat,n_cat,size(data,2));

for comp = 1:length(comparison)
    dsm(comparison(comp,1),comparison(comp,2),:) = data(comp,:); %3rd dim is time
    dsm(comparison(comp,2),comparison(comp,1),:) = data(comp,:);
end

if nargin>2 && nargout>1
    ii=1;
    avg_dsm = zeros(n_cat,n_cat,size(data,2)-(factor*2));
    
    for tp = 1+factor : size(dsm,3) - (factor)
        avg_dsm(:,:,ii) = mean(dsm(:,:,tp-factor:tp+factor),3);
        ii=ii+1;
    end
end
end