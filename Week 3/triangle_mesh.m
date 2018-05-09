rows=10;
cols=20;
tris=[];
for ii=1:rows-1
    i=ii-1;
    if(mod(ii,2)==1)
    for j =1:cols-1
        if(mod(j,2)==1)
            tris=[tris;cols*i+j cols*i+j+1 cols*(i+1)+j+1];
            tris=[tris;cols*i+j cols*(i+1)+j cols*(i+1)+j+1];
        else
            tris=[tris;cols*i+j cols*i+j+1 cols*(i+1)+j];
            tris=[tris;cols*(i+1)+j cols*(i)+j+1 cols*(i+1)+j+1];
        end
    end
    else 
    for j =1:cols-1
        if(mod(j,2)==0)
            tris=[tris;cols*i+j cols*i+j+1 cols*(i+1)+j+1];
            tris=[tris;cols*i+j cols*(i+1)+j cols*(i+1)+j+1];
        else
            tris=[tris;cols*i+j cols*i+j+1 cols*(i+1)+j];
            tris=[tris;cols*(i+1)+j cols*(i)+j+1 cols*(i+1)+j+1];
        end
    end
    end
end
[GX,GY] = meshgrid(1:rows,1:cols);
X=reshape(GX,[1,rows*cols]);
Y=reshape(GY,[1,rows*cols]);
triplot(tris,X,Y);