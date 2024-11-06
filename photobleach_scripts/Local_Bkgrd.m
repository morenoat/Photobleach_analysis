function [bkgImg,stdImg]=Local_Bkgrd(image,halfBox,Blocksearch,Border)


img=image;
Stepper=halfBox;
Blocksize=Blocksearch;
Rowcent=Border:Stepper:512-Border;
Colcent=Border:Stepper:256-Border;
bkgimg=zeros(512,256);
bkgimgSTD=zeros(512,256);

%tracker=zeros(512,256);
%C=1;
for jj=1:2:size(Rowcent,2)
    for ff=1:2:size(Colcent,2)
        [BXxlow,BXxhi,BXylow,BXyhi]=AOI_Limits([Rowcent(1,jj) Colcent(1,ff)],Blocksize/2);
        BXwndw=double(img(BXxlow:BXxhi,BXylow:BXyhi));
        A=BXwndw(:);
        zFactor = 1; % or whatever you want.
        
        outliers=1;
        while any(outliers)==1
            stdDev = std(A(:)); % Compute standard deviation
            medianValue = median(A(:)); % Compute median
            %Create a binary map of where outliers live
            outliers = abs(A-medianValue) > (zFactor * stdDev);
            A=A(outliers==0);
            stdDev = std(A(:)); % Compute standard deviation
            medianValue = median(A(:)); % Compute median
            %Create a binary map of where outliers live
            outliers = abs(A-medianValue) > (zFactor * stdDev);
        end
        [xlow,xhi,ylow,yhi]=AOI_Limits([Rowcent(1,jj) Colcent(1,ff)],Stepper);
        bkgimg(xlow:xhi,ylow:yhi)=medianValue;
        bkgimgSTD(xlow:xhi,ylow:yhi)=stdDev;
        %tracker(xlow:xhi,ylow:yhi)=C;
        %C=C+1;
        
    end
    
end

bkgImg=bkgimg;
stdImg=bkgimgSTD;

end

