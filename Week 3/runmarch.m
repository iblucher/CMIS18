clear all;
close all;
I = read_bw( 'EG_WEB_logo.jpg');  % Read black and white image from file
phi = bw2phi(I);
cols=24;
rows=24;
tris=[];
[GX,GY] = meshgrid(1:cols,1:rows);
[M, N] =size(phi);
X=reshape(GX,[1,cols*rows]);
Y=reshape(GY,[1,cols*rows]);
for ii=1:cols-1
    i=ii-1;
    if(mod(ii,2)==1)
    for j =1:rows-1
        if(mod(j,2)==1)
            tris=[tris;rows*i+j rows*i+j+1 rows*(i+1)+j+1];
            tris=[tris;rows*i+j rows*(i+1)+j rows*(i+1)+j+1];
        else
            tris=[tris;rows*i+j rows*i+j+1 rows*(i+1)+j];
            tris=[tris;rows*(i+1)+j rows*(i)+j+1 rows*(i+1)+j+1];
        end
    end
    else 
    for j =1:rows-1
        if(mod(j,2)==0)
            tris=[tris;rows*i+j rows*i+j+1 rows*(i+1)+j+1];
            tris=[tris;rows*i+j rows*(i+1)+j rows*(i+1)+j+1];
        else
            tris=[tris;rows*i+j rows*i+j+1 rows*(i+1)+j];
            tris=[tris;rows*(i+1)+j rows*(i)+j+1 rows*(i+1)+j+1];
        end
    end
    end
%    triplot(tris,X,Y);
%    pause(0.1);
end
%%
[GX,GY] = meshgrid(1:cols,1:rows);
[M, N] =size(phi);
X=reshape(GX,[1,cols*rows]);
X=X/cols*M;
Y=reshape(GY,[1,cols*rows]);
Y=Y/rows*N;
[NX,NY,newtris]=marching(Y,X,tris,phi,cols,rows);
%%
figure();
hold on
%imagesc(phi); colormap bone;
%triplot(tris,Y,X,"r");
triplot(newtris,NX,NY, 'red');
xlabel('x coordinates');
ylabel('y coordinates');
axis([0 165 0 105]);
hold off