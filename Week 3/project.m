function [X,Y] = project(X,Y,phi,GX,GY)
d=interp2(GX,GY,phi,X,Y);
[dx,dy]=gradient(phi);
nx=interp2(GX,GY,dx,X,Y);
ny=interp2(GX,GY,dy,X,Y);
dx=d.*nx;dy=d.*ny;
X(d>0)=X(d>0)-dx(d>0);
Y(d>0)=Y(d>0)-dy(d>0);
end

