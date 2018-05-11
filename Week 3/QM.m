function [shape,ssize]=QM(verts,tris)
size(verts)
shape=zeros([size(tris,1),1]);
ssize=zeros([size(tris,1),1]);
for i=1:size(tris,1)
p1=verts(tris(i,1),:);
p2=verts(tris(i,2),:);
p3=verts(tris(i,3),:);
l1=sqrt((p2(1)-p3(1))^2+(p2(2)-p3(2))^2);
l2=sqrt((p1(1)-p3(1))^2+(p1(2)-p3(2))^2);
l3=sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2);
x1=p1(1);x2=p2(1);x3=p3(1);
y1=p1(2);y2=p2(2);y3=p3(2);
shape(i)=(.5*((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1)))/(l1*l2*l3);
ssize(i)=1/max([l1,l2,l3]);
end