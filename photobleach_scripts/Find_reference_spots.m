function [Spots]=Find_reference_spots(Image,border,SpotSeparator,position_threshold,Bgcut,LargeBox,area,boxwidth)

        Iraw=Image;%image for spot indentification
        spotbx=SpotSeparator;%minimum pixels between spots
        Thresh=position_threshold;%maximum pixels between centroid center and 2D gaussian center
        halfBox=area/2;
        se = strel('disk',boxwidth,0);%structural element 
        MaskSpot = se.getnhood();
    

        [bkgImg,~]=local_background(Iraw,halfBox,LargeBox,border);%make background image by averaging local intensity
        bkgImg(bkgImg==0)=median(bkgImg(:));%replace background values = 0 with median background intensity 
        Isub=Iraw-bkgImg;%subtract local background 
        Ifit=relnoise(Isub,7,3,'disk');%local background subtraction
        [bestbind,~,~,~]=curateList_simple(Ifit(:),1,4);%remove outliers from intensity 
        bgThresh=prctile(bestbind,Bgcut);
        
        dpkRaw=double(peak_find_fast(Ifit,bgThresh));%FastPeakFind(d, thres, filt ,edg, res, fid)
        dpkRaw=[dpkRaw(1:2:end),dpkRaw(2:2:end)]; %xy coordniate of donor channel
        doutcent=cntrd(Ifit,round(dpkRaw),5);%calculates the centroid of bright spots to sub-pixel accuracy.
       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Remove Molecules that are too close to the border                  
        doutcent=doutcent(doutcent(:,2)>border+halfBox & doutcent(:,2)<512-(border+halfBox) & doutcent(:,1)>border+halfBox & doutcent(:,1)<256-(border+halfBox),:);% remove points that are within border
        [doutdat,dimg]=fitspot_2D_gaussian(Iraw,Ifit,bkgImg,doutcent(:,1:2),size(MaskSpot,1),MaskSpot);%Model spots as 2D gaussian for the shape and intensity of the spots
        doutdat(:,8)=0;%add selection boolean
        doutdat1=doutdat(doutdat(:,8)==0,:);
        dimg1=dimg;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        if size(doutdat,1)>1
            dDistmat=pdist2(doutdat(:,1:2),doutcent(:,1:2)); %calculate the pairwise distances between points in dchn 
            [dval,~]=sort(dDistmat,2); %sort the distances down columns 
            doutdat=doutdat(dval(:,2)>spotbx,:); %remove spots that are inside distance threshold leading to 'spurious detection'
            dimg=dimg(:,:,dval(:,2)>spotbx); %remove spots that are inside distance threshold leading to 'spurious detection'
        end
        
        if size(doutdat,1)>1
            dDistmat=pdist2(doutdat(:,1:2),doutcent(:,1:2)); %calculate the pairwise distances between points in dchn 
            [ddval,~]=sort(dDistmat,2); %sort the distances down columns 
            doutdat=doutdat(ddval(:,1)<Thresh,:); %remove spots where center drifts more than 1.25 pixels from centriod
            dimg=dimg(:,:,ddval(:,1)<Thresh); %remove spots that are inside distance threshold leading to 'spurious detection'
        end

        Spots.Centroid=doutcent;%Spots centers Identified using Centroid
        Spots.Psf=doutdat1;%2D guassian centers identified 
        Spots.Curate=doutdat;%curated spot list removing spots too close and dirfting too far between 2D and centroid
        Spots.PsfImage=dimg1;%image of 2D fit spots
        Spots.CurateImage=dimg;%image of curated spots
        Spots.BackGroundImg=bkgImg;
        Spots.SubImg=Isub;
        Spots.fitImg=Ifit;
end