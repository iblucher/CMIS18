function phi = bw2phi(I)
% BW2PHI  Convert black and white image to a signed distance field
% input:
%   I - The black and white image
% output
%   phi - The signed distance field 
% Copyright 2009, Kenny Erleben, DIKU.

phi = bwdist(I) - bwdist(~I);
ind = find(phi>0);
phi(ind) = phi(ind)-0.5;
ind = find(phi<0);
phi(ind) = phi(ind)+0.5;

end