function [SpotFit,dBndpss]=integrate_cy3_molecules(Image,Spot_local,sen)

        Iraw=Image;
        IsubB=relnoise(Iraw,7,3,'disk');%local background subtraction
        I_bpass = bpass(IsubB,1,5);
        Spot_local_B=Spot_local;
        Spot_local_B(:,3)=1;

%%%% Integrate spot intensity and save image of spots
    ElementStructure=sen;
    pts=Spot_local_B;
    tmp_array=zeros(size(pts,1),5);
    tnyBxsize=size(ElementStructure,1); % structure element size
    tnyEl=ElementStructure; %structuring element%
    bkgEl=abs(tnyEl-1);%
    bandpassimg=zeros(tnyBxsize,tnyBxsize,size(pts,1));%array to store spot images

    parfor ff=1:size(pts,1)
        pc2=zeros(5,1);
        aoix=pts(ff,2);
        aoiy=pts(ff,1);
        pc2(1,1)=pts(ff,1);
        pc2(2,1)=pts(ff,2);
        [xlowb,xhib,ylowb,yhib]=AOI_Limits([aoix aoiy],tnyBxsize/2);% define area to integrate
        wndftb=Iraw(xlowb:xhib,ylowb:yhib);%"cut out" spot area to integrate
        wndbp=I_bpass(xlowb:xhib,ylowb:yhib);%"cut out" bandpass spot area to visualize
        boxb = wndftb.*tnyEl;%select integration based on structural element 
        ptIb = sum(boxb(:));% integrate spot intensity
        Bkg = (wndftb.*bkgEl);% select background area
        Bkg = median(Bkg(Bkg(:)~=0));% determine median background
        pc2(3,1)=ptIb-(Bkg*sum(tnyEl(:)));% background subtracted spot intensity
        pc2(4,1)=0;%spot size not calculated
        pc2(5,1)=0;%spot distance not calculated
        tmp_array(ff,:)=pc2'; 
        bandpassimg(:,:,ff)=wndbp; %save bkgrd image of point for comparison
    end
    SpotFit=tmp_array;
    dBndpss=bandpassimg;

 
end



