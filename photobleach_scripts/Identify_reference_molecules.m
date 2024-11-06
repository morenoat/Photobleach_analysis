function []=Identify_reference_molecules(ngimages, exclude, outf,boxsize,DAdistance) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%

mkdir(pwd,outf);
dataout=[pwd '/' outf];
numngs=ngimages;
border=exclude;
a=pwd;
HighCutoff=99.5;%set high cutoff for removing outliers from distribution
LowCutoff=1;%set low cutoff for removing outliers from distribution
boxwidth = 4; % Width of box to draw around points
area=boxsize;%area of box to define the spot 12 by 12 pixels
prdist=DAdistance;% maximum distance between pairs of spots in Cy3 and Cy5 channel 
disThres=5;%separation between spots
Guass_Drift=3;%seperation between centroid and PSF
BGcut=85;%minimum level above babkground for spot detect
BoxDefineLocalBkgrd=40;%define box to determine local background in an image
NumSizeDeviate=5;% fold of StdDev to keep
cy3_normalization=50000;%scale Cy3 integrated intensity of the spots
cy5_normalization=50000;%scale Cy5 integrated intensity of the spots








%% Initialize Structres To Store Data
% Create structure array to store data from the different nanogird images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

nanogrids = arrayfun(@(K) struct('RawDch',zeros(512,216),'RawAch',zeros(512,216),'apk',[],'dpk',[],...
    'aPsf',[],'dPsf',[],'aPsfcrd',[],'dPsfcrd',[],'border',border), 1:numngs, 'UniformOutput',0);
      
    %%%% Split FOV into donor and Accpetor channels%%%%%%%
    for imind=1:numngs
        im=double(imread([a,'\ng' num2str(imind) '.tif']));
        RawDch = im(:,1:256); % donor half of the image
        RawAch= im(:,257:end); % acceptor half of the image
        
        %dchRaw = RawDch./65535; % donor half of the image
        %achRaw= RawAch./65535; % acceptor half of the image
        
        
        nanogrids{imind}.RawDch=RawDch;
        nanogrids{imind}.RawAch=RawAch; 
        %nanogrids{imind}.dchRaw=dchRaw;
        %nanogrids{imind}.achRaw=achRaw; 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Identify Molecules based on intensity in D-Channel      %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %[dSpots]=Findspots_BG(dchRaw,border,7,1.75,40,area,intMask);
        [dSpots]=Find_spots(RawDch,border,disThres,Guass_Drift,BGcut,BoxDefineLocalBkgrd,area,boxwidth);

        nanogrids{imind}.dpk=dSpots.Centroid;
        nanogrids{imind}.dPSFoutPre=dSpots.Psf;
        nanogrids{imind}.doutbx=dSpots.Curate;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Identify Molecules based on intensity in A-Channel      %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [aSpots]=Find_spots(RawAch,border,disThres,Guass_Drift,BGcut,BoxDefineLocalBkgrd,area,boxwidth);
        
        nanogrids{imind}.apk=aSpots.Centroid;
        nanogrids{imind}.aPSFoutPre=aSpots.Psf;
        nanogrids{imind}.aoutbx=aSpots.Curate;
    
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Determine the Intensity and spot shape distribution     %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cy3_all_molecules=[];
    cy5_all_molecules=[];
    for jj=1:numngs
        cy3_all_molecules=vertcat(cy3_all_molecules,nanogrids{jj}.doutbx(:,:));
        cy5_all_molecules=vertcat(cy5_all_molecules,nanogrids{jj}.aoutbx(:,:));
    end 
    cy3_all_molecules(:,12)=cy3_all_molecules(:,11)./cy3_normalization;%scale intensity of Cy3 molecules
    cy5_all_molecules(:,12)=cy5_all_molecules(:,11)./cy5_normalization;%scale intensity of Cy5 molecules
 
    [dbestMol, dremoveMol,dI_mean,~,~,dIstd]=curateList(cy3_all_molecules,12,3,HighCutoff,LowCutoff);
    dImx=prctile(cy3_all_molecules(:,12),HighCutoff);%determine high percentile cutoff for dchan spot selection
    %dImn=dLowCut;%determine low percentile cutoff for dchan spot selection
    dImn=dI_mean-(dIstd*1.0);%determine low percentile cutoff for dchan spot selection
    [dbestMolb, dremoveMolb,dradXmn,~,~,dradXsd]=curateList(dbestMol,5,3,HighCutoff,LowCutoff);% determine median spot radius X direction
    dradXsd=dradXsd*NumSizeDeviate;
    max_dX_radius=dradXmn+(dradXsd);
    min_dX_radius=dradXmn-(dradXsd/2);
    dremoveMol=vertcat(dremoveMol,dremoveMolb);
    [dbestMolfin, dremoveMolc,dradYmn,~,~,dradYsd]=curateList(dbestMolb,6,3,HighCutoff,LowCutoff);% determine median spot radius Y direction
    dremoveMol=vertcat(dremoveMol,dremoveMolc);
    dradYsd=dradYsd*NumSizeDeviate;
     max_dY_radius=dradYmn+(dradYsd);
    min_dY_radius=dradYmn-(dradYsd/2);
    
    [abestMol, aremoveMol,aI_mean,~,~,aIstd]=curateList(cy5_all_molecules,12,3,HighCutoff,LowCutoff);
    aImx=prctile(cy5_all_molecules(:,12),HighCutoff);%determine high percentile cutoff for Cy5 spot selection
    %aImn=aLowCut;%determine low percentile cutoff for dchan spot selection
    aImn=aI_mean-(aIstd*1.5);%determine low percentile cutoff for dchan spot selection
    [abestMolb, aremoveMolb,aradXmn,~,~,aradXsd]=curateList(abestMol,5,3,HighCutoff,LowCutoff);% determine median spot radius X direction
    aremoveMol=vertcat(aremoveMol,aremoveMolb);
    aradXsd=aradXsd*NumSizeDeviate;
    max_aX_radius=aradXmn+(aradXsd);
    min_aX_radius=aradXmn-(aradXsd/2);


    [abestMolfin, aremoveMolc,aradYmn,~,~,aradYsd]=curateList(abestMolb,6,3,HighCutoff,LowCutoff);% determine median spot radius Y direction
    aremoveMol=vertcat(aremoveMol,aremoveMolc);
    aradYsd=aradYsd*NumSizeDeviate;
    
    max_aY_radius=aradYmn+(aradYsd);
    min_aY_radius=aradYmn-(aradYsd/2);
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%            Plot selection information                   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %plotIvsRad_color_select(Intensity,Radius,max_intensity,min_Intensity,max_radius,min_radius,Xaxis,Yaxis,color_all_points,color_select_points,plotName)
    
    
     plotIvsRad_color_select(cy3_all_molecules(:,12),cy3_all_molecules(:,5),dImx,dImn,max_dX_radius,min_dX_radius,'Intensity',...
         'X Radius (Pixels)',25,'b',[0.8500 0.3250 0.0980],'Cy3 Intensity vs Xradius',dataout,outf,1);
    plotIvsRad_color_select(cy3_all_molecules(:,12),cy3_all_molecules(:,6),dImx,dImn,max_dY_radius,min_dY_radius,'Intensity',...
        'Y Radius (Pixels)',25,'b',[0.8500 0.3250 0.0980],'Cy3 Intensity vs Yradius',dataout,outf,1);
    plotIvsRad_color_select(cy5_all_molecules(:,12),cy5_all_molecules(:,5),aImx,aImn,max_aX_radius,min_aX_radius,'Intensity',...
        'X Radius (Pixels)',25,'b',[0.8500 0.3250 0.0980],'Cy5 Intensity vs Xradius',dataout,outf,1);
    plotIvsRad_color_select(cy5_all_molecules(:,12),cy5_all_molecules(:,6),aImx,aImn,max_aY_radius,min_aY_radius,'Intensity',...
        'Y Radius (Pixels)',25,'b',[0.8500 0.3250 0.0980],'Cy5 Intensity vs Yradius',dataout,outf,1);
     
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%            Save selection information                   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
    nanogrids{numngs+1}.dPsfTot=cy3_all_molecules;
    nanogrids{numngs+1}.aPsfTot=cy5_all_molecules;
    nanogrids{numngs+1}.dIntMaxcut=dImx;
    nanogrids{numngs+1}.dIntMincut=dImn;
    nanogrids{numngs+1}.d_mean_intensity=dI_mean;
    nanogrids{numngs+1}.dmeanXrad=dradXmn;
    nanogrids{numngs+1}.max_dX_radius=max_dX_radius;
    nanogrids{numngs+1}.min_dX_radius=min_dX_radius;
    nanogrids{numngs+1}.max_dY_radius=max_dY_radius;
    nanogrids{numngs+1}.min_dY_radius=min_dY_radius;    
    
    
    nanogrids{numngs+1}.dstdXrad=dradXsd;
    nanogrids{numngs+1}.dmeanYrad=dradYmn;
    nanogrids{numngs+1}.dstdYrad=dradYsd;
    nanogrids{numngs+1}.dbestMolfin=dbestMolfin;
    nanogrids{numngs+1}.dremoveMol=dremoveMol;
    
 
    nanogrids{numngs+1}.aIntMaxcut=aImx;
    nanogrids{numngs+1}.aIntMincut=aImn;
    nanogrids{numngs+1}.a_mean_intensity=aI_mean;

    nanogrids{numngs+1}.max_aX_radius=max_aX_radius;
    nanogrids{numngs+1}.min_aX_radius=min_aX_radius;
    nanogrids{numngs+1}.max_aY_radius=max_aY_radius;
    nanogrids{numngs+1}.min_aY_radius=min_aY_radius;  

    nanogrids{numngs+1}.ameanXrad=aradXmn;
    nanogrids{numngs+1}.astdXrad=aradXsd;
    nanogrids{numngs+1}.ameanYrad=aradYmn;
    nanogrids{numngs+1}.astdYrad=aradYsd;
    nanogrids{numngs+1}.abestMolfin=abestMolfin;
    nanogrids{numngs+1}.aremoveMol=aremoveMol; 
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     Spot selection based on intensity & Radius          %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    
    for xx=1:numngs
        
        dDat=nanogrids{xx}.doutbx;
        %dImg=nanogrids{xx}.dimgbx;
        aDat=nanogrids{xx}.aoutbx;
       % aImg=nanogrids{xx}.aimgbx
        dDat(:,12)=dDat(:,11)./cy3_normalization;
        aDat(:,12)=aDat(:,11)./cy5_normalization;
        dDatB=dDat(dDat(:,12)>dImn & dDat(:,12)<dImx & dDat(:,5)>=min_dX_radius & dDat(:,5)<=max_dX_radius ...
            & dDat(:,6)>=min_dY_radius & dDat(:,6)<=max_dY_radius,:);
            %keep points with intensity and radius within acceptable range
       % dImgB=dImg(:,:,dDat(:,3)>dImn & dDat(:,3)<dImx & dDat(:,5)>=dradXmn-dradXsd & dDat(:,5)<=dradXmn+dradXsd ...
       % & dDat(:,6)>=dradYmn-dradYsd & dDat(:,6)<=dradYmn+dradYsd);
        %keep points with intensity and radius within acceptable range
        aDatB=aDat(aDat(:,12)>aImn & aDat(:,12)<aImx & aDat(:,5)>=min_aX_radius & aDat(:,5)<=max_aX_radius ...
             & aDat(:,6)>=min_aY_radius & aDat(:,6)<=max_aY_radius,:);
        %keep points with intensity and radius within acceptable range
        %aImgB=aImg(:,:,aDat(:,3)>aImn & aDat(:,3)<aImx & aDat(:,5)>=aradXmn-aradXsd & aDat(:,5)<=aradXmn+aradXsd ...
            % & aDat(:,6)>=aradYmn-aradYsd & aDat(:,6)<=aradYmn+aradYsd);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%    Remove markers that are not present in both channels     %%%
    %%%                                                             %%%  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
        if ~isempty(dDatB) && ~isempty(aDatB)
            
            dchDistmat=pdist2(dDatB(:,1:2),aDatB(:,1:2)); %calculate the pairwise distances between points in dchn and achn
            [ddval,~]=sort(dchDistmat,2);%sort the distances down columns 
            [aaval,~]=sort(dchDistmat,1);
            fnDch=dDatB(ddval(:,1)<prdist,:); % Save points that are paired
            %fnDimg=dImgB(:,:,ddval(:,1)<prdist);
            fnAch=aDatB(aaval(1,:)'<prdist,:);% Save points that are paired
            %fnAimg=aImgB(:,:,aaval(1,:)'<prdist);
        
        end
        if isempty(dDatB) && isempty(aDatB)
            fnAch=[];
            fnDch=[];
        end
        
        if ~isempty(fnDch) && ~isempty(fnAch)

            chnDistmat=pdist2(fnDch(:,1:2),fnAch(:,1:2)); %calculate the pairwise distances between points in dchn and achn
            [~,dixA]=sort(chnDistmat,2);
            achCrd=fnAch(dixA(:,1),:);
            %aimgPst=fnAimg(:,:,dixA(:,1));
            valDist=pdist2(achCrd(:,1:2),fnDch(:,1:2));
            [~,dixD]=sort(valDist,2);
            dchCrd=fnDch(dixD(:,1),:);
            %dimgPst=fnDimg(:,:,dixD(:,1));
            nanogrids{xx}.achCrd=achCrd;
            nanogrids{xx}.dchCrd=dchCrd;
            %nanogrids{xx}.aimgPst=aimgPst;
           % nanogrids{xx}.dimgPst=dimgPst;
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%                                           %%%
        %%%              Save "Bad Points"            %%%  
        %%%                                           %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dpreDat=nanogrids{xx}.dPSFoutPre;
        apreDat=nanogrids{xx}.aPSFoutPre;
        
        dPre=round(dpreDat(:,1:2));
        dPST=round(dchCrd(:,1:2));
        Lid = ismember(dPre,dPST,'rows');
        baddpt=dpreDat(Lid==0,:);
        dpreDat(Lid==0,9)=0;
        dpreDat(Lid==1,9)=1;
        aPre=apreDat(:,1:2);
        aPST=achCrd(:,1:2);
        Lia = ismember(aPre,aPST,'rows');
        badapt=apreDat(Lia==0,:);
        apreDat(Lia==0,9)=0;
        apreDat(Lia==1,9)=1;        

        nanogrids{xx}.badDpt=baddpt;
        nanogrids{xx}.badApt=badapt;
        nanogrids{xx}.dPSFoutPre=dpreDat;
        nanogrids{xx}.aPSFoutPre=apreDat;
        nanogrids{xx}.post_curate_Cy3=dDatB;
        nanogrids{xx}.post_curate_Cy5=aDatB;


        else
            achCrd=[];
            dchCrd=[];
            aaval=[];
            ddval=[];
            nanogrids{xx}.achCrd=achCrd;
            nanogrids{xx}.dchCrd=dchCrd;
            nanogrids{xx}.aval=aaval;
            nanogrids{xx}.dval=ddval;
            nanogrids{xx}.badDpt=[];
            nanogrids{xx}.badApt=[];
        end
    end
    %close all      
    %clearvars -except nanogrids a numngs border imind outf boxwidth dataout
                   

%Combine all marker points from different nanogrid images

    totDpt=[];
    totApt=[];
    for x=1:numngs
        totDpt=vertcat(totDpt,nanogrids{x}.dchCrd); 
        totApt=vertcat(totApt,nanogrids{x}.achCrd);
    end
    nanogrids{numngs+1}.totDpt=totDpt;
    nanogrids{numngs+1}.totApt=totApt;
    
save(fullfile(dataout,outf),'nanogrids','-v7.3','-nocompression' );
end 
    

    

    
    
    
    
    
    
    
    
    