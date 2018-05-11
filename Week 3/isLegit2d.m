function [res] = isLegit2d(n,TR,i,phi)
%to validate whether a triangle is to be kept, check wheter n random points
%inside triangle ps is > in phi
count=0;

[M, N] =size(phi);
[GX,GY] = meshgrid(1:N,1:M);
res=1;
while count < n
    w=rand(1,3);
    w=w./sum(w);
    p=barycentricToCartesian(TR,i,w);
    if 0>interp2(GX,GY,phi,p(1),p(2))
       count=count+1;
    else
       res=0;
       break
    end
end

end

