function [spot_intensity_position_size_array,spot_image_output]=integrate_cy5_molecules(Image,Spot_local,border,Tehter_Limit,Bgcut,LargeBox,area,sen)
%This function will ID spots in the cy5 channel based on maximum intensity.
%It will then keep spot location that are within a distance threshold of
%cy3 molecules. The function will cut out window from image (fitimage)
%defined by box_size. For each spot, the background integrated intensity is calculated and the spot is fit to
% a 2D gaussian
%The output (outdat) is an array each row represents a spot columns are:
%[x-position; y-position; Amplitude; background;sigma]
        Iraw=Image;
        Thresh=Tehter_Limit;
        halfBox=area/2;
        [bkgImg,~]=Local_Bkgrd(Iraw,halfBox,LargeBox,border);
        bkgImg(bkgImg==0)=median(bkgImg(:));
        Isub=Iraw-bkgImg;%subtract local background 
        Ifit=relnoise(Isub,7,3,'disk');%local background subtraction
        IsubB=relnoise(Iraw,7,3,'disk');%local background subtraction
        I_bpass = bpass(IsubB,1,5);
        Allval=Ifit(:);
        [bestbind,~,~,~]=curateList_simple(Allval,1,4);
        bgThresh=prctile(bestbind,Bgcut);
        dpkRaw=double(peak_find_fast(Ifit,bgThresh));%FastPeakFind(d, thres, filt ,edg, res, fid)
        dpkRaw=[dpkRaw(1:2:end),dpkRaw(2:2:end)]; %xy coordniate of donor channel
        doutcent=cntrd(Ifit,round(dpkRaw),5);%calculates the centroid of bright spots to sub-pixel accuracy.
        Spot_local_B=Spot_local;
        Spot_local_B(:,3)=0;
        if size(doutcent,1)>1
            dDistmat=pdist2(doutcent(:,1:2),Spot_local(:,1:2)); %calculate the pairwise distances between points in dchn
            [dval,didx]=sort(dDistmat,2); %sort the distances down columns
            Spot_local_B(didx(dval(:,1)<Thresh),3)=1;
            Spot_local_B(didx(dval(:,1)<Thresh),1:2)=doutcent(dval(:,1)<Thresh,1:2);
        end
       

pts=Spot_local_B;
ptsA=Spot_local;
bx=area;
tmp_data_array=zeros(size(pts,1),5);

tnyBx=size(sen,1); %size of box for structure element
tnyEl=sen; %structuring element
bkgEl=abs(tnyEl-1);% background structure element
spot_bandpass_images=zeros(tnyBx,tnyBx,size(pts,1));

parfor ff=1:size(pts,1)
    pc=[];
    tmppt=ptsA(ff,:);
    aoix=pts(ff,2);
    aoiy=pts(ff,1);
    [xlow,xhi,ylow,yhi]=AOI_Limits([aoix aoiy],bx/2);
    wndft=Iraw(xlow:xhi,ylow:yhi);
    %wndbp=BandPass(xlow:xhi,ylow:yhi);
    mxI=max(wndft(:));
    medI=mean(wndft(:));

    % Full usage: [ params, exitflag ] = psfFit_Image( img, param_init, param_optimizeMask, useIntegratedGauss, useMLErefine, hWinSize, global_init )
    %dPsf=psfFit_Image(dchAve, dpks', [1,1,1,1,1,1,1]);%Fit peaks to 2D gaussian for greater precission [xpos,ypos,A,BG,sigma_x,sigma_y,angle]
%The output (pc2) is an array each row represents a spot columns are:
%[x-position; y-position; Amplitude; background;sigma]
    pc2=psfFit_Image(wndft, [bx/2;bx/2;mxI-medI;medI;bx/4], [1,1,1,1,1,0,0]);
    pc(1,1)=ylow+pc2(1,1)-1;
    pc(2,1)=xlow+pc2(2,1)-1;

    %%%Test levels using smaller structuring element
    [xlowb,xhib,ylowb,yhib]=AOI_Limits([aoix aoiy],tnyBx/2);
    wndftb=Iraw(xlowb:xhib,ylowb:yhib);%%%%%%%%%%%%%%%%%%%%
    wndbp=I_bpass(xlowb:xhib,ylowb:yhib);
    boxb = wndftb.*tnyEl;
    ptIb = sum(boxb(:));%%%%%%%%%%%%%%%%%%%%%%%
    Bkg = (wndftb.*bkgEl);%%%%%%%%%%%%%%%%%%%%%%%
    Bkg = median(Bkg(Bkg(:)~=0));%%%%%%%%%%%%%%%%
    pc(3,1)=ptIb-(Bkg*sum(tnyEl(:)));%Spot Intensity
    pc(4,1)=pc2(5,1);%pc3(5,1)=spot size
    pc(5,1)=sqrt((pc(1,1)-tmppt(1,1))^2 +(pc(2,1)-tmppt(1,2))^2);%Distance from tether = sqrt(fit(1,1)-(datapt_original(1,1))^2+fit(2,1)-(datapt_original(1,2))^2);
    tmp_data_array(ff,:)=pc';
    spot_bandpass_images(:,:,ff)=wndbp; %save bkgrd image of point for comparison
end
spot_intensity_position_size_array=tmp_data_array;
spot_image_output=spot_bandpass_images;





end