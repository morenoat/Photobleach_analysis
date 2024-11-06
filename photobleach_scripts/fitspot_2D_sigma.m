function [outdat,pointimg]=fitspot_2D_sigma(fitimage,subimg,bkgimg,datapt,box_size,SpotMask)
%Function 2Dfit(fitimage,datapt,box_size)
%This funtion will cut out window from image (fitimage) defined by box_size and model
%the intensity as a 2D gaussian for all xy cooridanates given by datapt
%The output (outdat) is an array each row represents a spot columns are:
%[x-position; y-position; Amplitude; background; x-sigma; y-sigma; angle; convergence ] 
    img=subimg;
    imgbp=bkgimg;
    pts=datapt;
    bx=box_size;
    tmp=zeros(size(pts,1),11);
    tmpimg=zeros(bx,bx,size(pts,1));
    fitimg=fitimage;
    sen=SpotMask;
    IntNum=sum(sen(:));
    BGbx=sen-1;
    BGbx(BGbx==-1)=1;
    
    for ff=1:size(pts,1)
        pc2=zeros(11,1);
        aoix=pts(ff,2);
        aoiy=pts(ff,1);
        [xlow,xhi,ylow,yhi]=AOI_Limits([aoix aoiy],bx/2);
        wndw=img(xlow:xhi,ylow:yhi);
        wndwbg=imgbp(xlow:xhi,ylow:yhi);
        wndft=fitimg(xlow:xhi,ylow:yhi);
        box = wndft.*sen;
        boxBg=wndft.*BGbx;
        ptI = sum(box(:));
        bgI = median((boxBg(boxBg(:)~=0)));%median background intensity
% Full usage: [ params, exitflag ] = psfFit_Image( img, param_init, param_optimizeMask, useIntegratedGauss, useMLErefine, hWinSize, global_init )
%dPsf=psfFit_Image(dchAve, dpks', [1,1,1,1,1,1,1]);%Fit peaks to 2D gaussian for greater precission [xpos,ypos,A,BG,sigma_x,sigma_y,angle] 
        pc2=psfFit_Image(wndft, [bx/2;bx/2;ptI;bgI;bx/4], [1,1,1,1,1,0,0]);
        pc2(6,1)=pc2(1,1);
        pc2(7,1)=pc2(2,1);
        pc2(1,1)=ylow+pc2(1,1)-1;
        pc2(2,1)=xlow+pc2(2,1)-1;
        pc2(9,1)=ptI;
        pc2(10,1)=bgI;
        pc2(11,1)=bgI*(IntNum);
        pc2(8,1)=ptI-(bgI*(IntNum));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        tmp(ff,:)=pc2'; 
        tmpimg(:,:,ff)=wndw; %save image of point for comparison
    end
    outdat=tmp;
    pointimg=tmpimg;
end

