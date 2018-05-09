function [NX,NY,newtris] = marching(X,Y,tris,phi)
newtris=[];
NX=X;
NY=Y;
[~,N]=size(NX);
tricount=size(tris);
tricount=tricount(1);
[M, NN] =size(phi);
[GX,GY] = meshgrid(1:NN,1:M);
for i=1:tricount
%    size(newtris)
    in1=0>interp2(GX,GY,phi,X(tris(i,1)),Y(tris(i,1)));
    in2=0>interp2(GX,GY,phi,X(tris(i,2)),Y(tris(i,2)));
    in3=0>interp2(GX,GY,phi,X(tris(i,3)),Y(tris(i,3)));
    inside=in1+in2+in3;
    if(inside==3)
        newtris=[newtris;tris(i,:)];
    elseif(inside==2)
        if(in1~=1)
            pout=tris(i,1);
            pin1=tris(i,2);
            pin2=tris(i,3);
        elseif(in2~=1)
            pout=tris(i,2);
            pin1=tris(i,3);
            pin2=tris(i,1);
        else
            pout=tris(i,3);
            pin1=tris(i,1);
            pin2=tris(i,2);
        end
            new1x=(X(pout)+X(pin1))/2;
            new1y=(Y(pout)+Y(pin1))/2;
            new2x=(X(pout)+X(pin2))/2;
            new2y=(Y(pout)+Y(pin2))/2;
            NX=[NX new1x new2x];
            NY=[NY new1y new2y];
            N=N+2;
            index1=N;
            index2=N-1;
            newtris=[newtris;pin1 pin2 index1];
            newtris=[newtris;pin1 index1 index2];
    elseif(inside==1)
        if(in1==1)
            pin=tris(i,1);
            pout1=tris(i,2);
            pout2=tris(i,3);
        elseif(in2==1)
            pin=tris(i,2);
            pout1=tris(i,3);
            pout2=tris(i,1);
        else
            pin=tris(i,3);
            pout1=tris(i,1);
            pout2=tris(i,2);
        end
            new1x=(X(pout1)+X(pin))/2;
            new1y=(Y(pout1)+Y(pin))/2;
            new2x=(X(pout2)+X(pin))/2;
            new2y=(Y(pout2)+Y(pin))/2;
            NX=[NX new1x new2x];
            NY=[NY new1y new2y];
            N=N+2;
            index1=N-1;
            index2=N;
            newtris=[newtris;pin index1 index2];
    else
        %no points inside dont care
    end
end

