function I = read_bw( filename )
% READ_BW - Read black and white image from a filename
% input: 
%  filename  - The filename of the image file to read
% output:
%  I         - The resulting black and white image
% Copyright 2009, Kenny Erleben, DIKU.
I = double( ~im2bw( imread( filename ) ) );
end
