function [legitTris] = delMesh(X,Y,phi)
T=delaunay(X,Y);
%triplot(T,X,Y);
hold off
legitTris=[];
tris=size(T);
TR = triangulation(T,[X;Y]');
for i=1:tris(1)
    p1=T(i,1);
    p2=T(i,2);
    p3=T(i,3);
    if(isLegit2d(10,TR,i,phi))
        legitTris=[legitTris;p1 p2 p3];
    end
end
end

