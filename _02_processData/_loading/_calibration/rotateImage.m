%{
------------HEADER-----------------
Objective          ::  Define coordinates in any domain.

INPUT VARS
im_origin          :: the original image
angleInDegree  :: The angle the image must be rotated

OUTPUT VARS
imagerot             :: the rotated image

Created by         :: O. Beltramo-Martin - ONERA/LAM
Creation date      :: 11/05/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function imagerot = rotateImage(im_origin,angleInDegree)
            switch mod(angleInDegree, 360)
                % Special cases
                case 0
                    imagerot = im_origin;
                case 90
                    imagerot = rot90(im_origin);
                case 180
                    imagerot = im_origin(end:-1:1, end:-1:1);
                case 270
                    imagerot = rot90(im_origin(end:-1:1, end:-1:1));                    
                otherwise
                    % Zero-pad the image
                    [Rows, Cols] = size(im_origin);
                    Diagonal = hypot(Rows,Cols);
                    RowPad = ceil(Diagonal - Rows) + 2;
                    ColPad = ceil(Diagonal - Cols) + 2;
                    im_pad = zeros(Rows+RowPad, Cols+ColPad);
                    im_pad(ceil(RowPad/2):(ceil(RowPad/2)+Rows-1),ceil(ColPad/2):(ceil(ColPad/2)+Cols-1)) = im_origin;

                    % Convert to radians and create transformation matrix                    
                    imagerot=zeros(size(im_pad)); % midx and midy same for both
                    angleInRad = angleInDegree*pi/180;
                    midx=ceil((size(im_pad,1)+1)/2);
                    midy=ceil((size(im_pad,2)+1)/2);
                    
                    for i=1:size(im_pad,1)
                        for j=1:size(im_pad,2)
                            
                            x= (i-midx)*cos(angleInRad)+(j-midy)*sin(angleInRad);
                            y=-(i-midx)*sin(angleInRad)+(j-midy)*cos(angleInRad);
                            x=round(x)+midx;
                            y=round(y)+midy;
                            
                            if (x>=1 && y>=1 && x<=size(im_pad,2) && y<=size(im_pad,1))
                                imagerot(i,j)=im_pad(x,y); % k degrees rotated image
                            end                            
                        end
                    end
                    imagerot = tools.crop(imagerot,size(im_origin));
            end
        end