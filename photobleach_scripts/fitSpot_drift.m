function [outdat,mappt]=fitSpot_drift(fitimage,subimg,bkgimg,datapt,area,SpotMask,sigthres)
%Function 2Dfit(fitimage,datapt,box_size)
%This funtion will cut out window from image (fitimage) defined by box_size and model
%the intensity as a 2D gaussian for all xy cooridanates given by datapt
%The output (outdat) is an array each row represents a spot columns are:
%[x-position; y-position; Amplitude; background; x-sigma; y-sigma; angle; convergence ] 
    img=subimg;
    imgbp=bkgimg;
    pts=datapt;
    bx=size(SpotMask,1);
    bxb=area;
    tmp=zeros(size(pts,1),11);
    fitimg=fitimage;
    sen=SpotMask;

    for ff=1:size(pts,1)
        pc2=zeros(11,1);
        aoix=pts(ff,2);
        aoiy=pts(ff,1);
        [xlow,xhi,ylow,yhi]=AOI_Limits([aoix aoiy],bx/2);
        [xlowb,xhib,ylowb,yhib]=AOI_Limits([aoix aoiy],bxb/2);
        wndwb=img(xlowb:xhib,ylowb:yhib);
        wndwbg=imgbp(xlow:xhi,ylow:yhi);
        wndft=fitimg(xlow:xhi,ylow:yhi);
        box = wndft.*sen;
        boxBg=wndwbg.*sen;
        ptI = sum(box(:));
        bgI = sum(boxBg(:));
        ptFit=[bx/2,bx/2];
        ptout=cntrd(wndwb,round(ptFit),3);%calculates the centroid of bright spots to sub-pixel accuracy.
        pc2(1,1)=ylow+ptout(1,1)-1;
        pc2(2,1)=xlow+ptout(1,2)-1;
        pc2(9,1)=ptI;
        pc2(10,1)=bgI;
        pc2(11,1)=ptI./bgI;
        
        if pc2(11,1)>sigthres && ylow+pc2(1,1)-1>bx+1 && ylow+pc2(1,1)-1<256-bx-1 && xlow+pc2(2,1)-1>bx+1 && xlow+pc2(2,1)-1<512-bx-1
            sinpt(1,1)=ylow+pc2(1,1)-1;
            sinpt(2,1)=xlow+pc2(2,1)-1;
        else
            sinpt(1,1)=pts(ff,1);
            sinpt(2,1)=pts(ff,2);
        end
  
        tmp(ff,:)=pc2'; 
        tmp(ff,:)=pc2';
        tmppt(ff,:)=sinpt';
    end
    outdat=tmp;
    mappt=tmppt;
end

