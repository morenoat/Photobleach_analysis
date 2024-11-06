function [segments, segtraj] = spot_selection_integration_frame_by_frame(imstack,infile,outf,fitgrid,varargin)

mkdir(pwd,[outf]);
dataout=[pwd '/' outf];
segtraj = {};
nanogrids=importdata([fitgrid '.mat']);
if ~exist('bfGetReader')
    % IF NECESSARY, CHANGE THIS TO THE CORRECT PATH.
    %addpath('C:\Program Files\MATLAB\R2018a\toolbox\bfmatlab');
end

% default correlation threshold with respect to maximum and minimum frame
% correlation values:
reader = bfGetReader(imstack);
% PARSE OPTIONAL INPUT ARGUMENTS
if ~isempty(varargin)
    argnum = 1;
    while argnum <= numel(varargin)
        switch varargin{argnum}
            case('segments')
                segments = varargin{argnum+1};
                argnum = argnum + 2;
            case('segmentfile')
                segmentfile = varargin{argnum+1};
                segments = open(segmentfile);
                segments = segments.segments;
                argnum = argnum + 2;
            case('MovieTime')
                timeMovie = varargin{argnum+1};
                argnum = argnum + 2;
        end
    end
end


inargs = importdata(infile);  % This will give trouble if your file is strictly numeric (it will import it as an array)
border = inargs.data(ismember(inargs.textdata, 'border'));
boxwidth = inargs.data(ismember(inargs.textdata, 'boxwidth'));
rectangle = inargs.data(ismember(inargs.textdata, 'rectangle'));
HighCutoff = inargs.data(ismember(inargs.textdata, 'HighCutoff'));
LowCutoff = inargs.data(ismember(inargs.textdata, 'LowCutoff'));
drkwindow = inargs.data(ismember(inargs.textdata, 'drkwindow'));
area = inargs.data(ismember(inargs.textdata, 'area'));
minimum_distance_between_centers = inargs.data(ismember(inargs.textdata, 'minimum_distance_between_centers'));
BGcut = inargs.data(ismember(inargs.textdata, 'BGcut'));
Gauss_Drift = inargs.data(ismember(inargs.textdata, 'Gauss_Drift'));
BoxDefineLocalBkgrd = inargs.data(ismember(inargs.textdata, 'BoxDefineLocalBkgrd'));
pixel_limit_from_cy3 = inargs.data(ismember(inargs.textdata, 'pixel_limit_from_cy3'));
fold_deviation_radius = inargs.data(ismember(inargs.textdata, 'fold_deviation_radius'));
number_Points_to_Fit = inargs.data(ismember(inargs.textdata, 'number_Points_to_Fit'));
cy3_scale = inargs.data(ismember(inargs.textdata, 'cy3_scale'));
number_of_frames_to_average = inargs.data(ismember(inargs.textdata, 'number_of_frames_to_average'));
max_distance_between_cy3_cy5 = inargs.data(ismember(inargs.textdata, 'number_of_frames_to_average'));
numFOVID=size(segments(:,1),1);
nseg=size(segments(:,1),1);
segtraj = {};
 
bar_wait = waitbar(0,'Loading Data','Name','ID fluorescent Molecules');




% mask that will be used to "cut out" spots later on for background
% determination:
se = strel('disk',boxwidth,0);
sen = se.getnhood();

close all;

for IDmol = 1:numFOVID  % Loop over segments in distinct fields of view

    IDstart=segments(IDmol,1);%Start point (Frame) to generate average image and ID tehter molecule
    MovieEnd=segments(IDmol,2);%End point (Frame) of trajectory


    %index of FOV start totTime is segments(:,1)
    %index of FOV stop totTime is segments(:,2)

    %time of FOV start in totTime is segments(:,3)
    %time of FOV stop  in totTime is segments(:,4)

    %index of FOV start in Cy3Time is segments(:,5)
    %index of FOV stop in Cy3Time is segments(:,6)

    %time of FOV start in Cy3Time is segments(:,7)
    %time of FOV stop  in Cy3Time is segments(:,8)

    %index of FOV start in Cy5Time is segments(:,9)
    %index of FOV stop in Cy5Time is segments(:,10)

    %time of FOV start in Cy5Time is segments(:,11)
    %time of FOV stop  in Cy5Time is segments(:,12)




    % Initialize data structure containing trajectories for this segment.
     segtraj{IDmol} = struct('segments',segments,'movieTime',timeMovie,'cy3_img_average',zeros(512,256),'cy5_img_average',zeros(512,256),...
       'overlaid_img_cy3_cy5',zeros(512,256,3),'id_cy3_spots',[],'id_cy5_spots',[],'MovieStart',[],'MovieEnd',[],...
       'SpotSelect',[],'offx',[],'offy',[],'drift_info',[],'cy3_final_spots',[],'cy5_final_spots',[],...
       'cy3_dark_final_roi',[],'cy5_dark_final_roi',[],'cy3_rejected_spots',[],'dark_roi_info',[],...
       'cy3_intensity',[],'cy3_dark_roi_intensity',[],'cy3_size',[],'cy3_dark_roi_size',[],...
       'cy3_distance',[],'cy3_dark_roi_distance',[],'cy5_intensity',[],'cy5_dark_roi_intensity',[],'cy5_size',[],'cy5_dark_roi_size',[],...
       'cy5_distance',[],'cy5_dark_roi_distance',[],'cy5_x_position',[],'cy5_y_position',[],...
       'cy5_dark_roi_x_position',[],'cy5_dark_roi_y_position',[]);
 

    cy3_img=zeros(512,256);
    cy5_img=zeros(512,256);
    a=0;
    b=0;
    gh=0;
    n=IDstart;
    while n<MovieEnd
        n=IDstart+gh;
        if timeMovie(n,2)==1 && a <= number_of_frames_to_average
            timg=double(bfGetPlane(reader,n));
            cy3_img = cy3_img + timg(:,1:256);
            a=a+1;
            gh=gh+1;
        elseif timeMovie(n,2)==2
            timg=double(bfGetPlane(reader,n));
            cy5_img = cy5_img + timg(:,257:512);
            b=b+1;
            gh=gh+1;
        else
            gh=gh+1;
        end

    end

 
    cy3_img_average = cy3_img./a;
    cy5_img_average = cy5_img./b;

    %% create rgb images to show overlap between between different fluorescent molecules
    noneArry=zeros(512,256);

    qc_Cy5=bpass(cy5_img_average,1,5);%bandpass filter Cy5 FOV
    Y = prctile(qc_Cy5(qc_Cy5(:)>0),70);%set Cy5 background make it the same across image
    qc_Cy5(qc_Cy5(:)<=Y)=Y;
    Yh=prctile(qc_Cy5(qc_Cy5(:)~=Y),97);%set Cy5 threshold for spot intensity
    qc_Cy5(qc_Cy5(:)>=Yh)=Yh;
    Cy5aa=rescale(qc_Cy5);%normalize the image
    qc_Cy3=cy3_img_average;
    YD=prctile(qc_Cy3(qc_Cy3(:)>0),70);%set Cy3 background make it the same across image
    qc_Cy3(qc_Cy3(:)<=YD)=YD;
    YDh=prctile(qc_Cy3(qc_Cy3(:)~=YD),97);%set Cy3 threshold for spot intensity
    qc_Cy3(qc_Cy3(:)>=YDh)=YDh;
    Cy3dd=rescale(qc_Cy3);%normalize the image


    %make rgb image with Cy5 FOV (red) and cy3 FOV (green)
    ImageV2(:,:,1)=Cy5aa;
    ImageV2(:,:,2)=Cy3dd.*.5;
    ImageV2(:,:,3)=noneArry;

    segtraj{IDmol}.cy3_img_average=cy3_img_average; %Save image for comparison
    segtraj{IDmol}.cy5_img_average=cy5_img_average; %Save image for comparison
    segtraj{IDmol}.overlaid_img_cy3_cy5=ImageV2; %Save image for comparison


   waitbar(.25,bar_wait,'ID molecules');




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%         Identify Molecules based on intensity           %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [id_cy3_spots]=Find_spots(cy3_img_average,border,minimum_distance_between_centers,Gauss_Drift,BGcut,BoxDefineLocalBkgrd,area,boxwidth);

    segtraj{IDmol}.id_cy3_spots=id_cy3_spots;

    [id_cy5_spots]=Find_spots(cy5_img_average,border,minimum_distance_between_centers,Gauss_Drift,BGcut,BoxDefineLocalBkgrd,area,boxwidth);
    segtraj{IDmol}.id_cy5_spots=id_cy5_spots;








end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Determine the Intensity and spot shape distribution     %%%
%%%             for Cy3 molecules                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cy3_all_molecules=[];
for jj=1:numFOVID
    cy3_all_molecules=vertcat(cy3_all_molecules,segtraj{jj}.id_cy3_spots.Curate);
end
cy3_all_molecules(:,12)=cy3_all_molecules(:,11)./cy3_scale;%scale intensity of Cy3 molecules
[dbestMol, ~,dI_mean,~,~,dIstd]=curateList(cy3_all_molecules,12,3,HighCutoff,LowCutoff);
cy3_max_intensity=prctile(cy3_all_molecules(:,12),HighCutoff);%determine high percentile cutoff for dchan spot selection
cy3_min_intensity=prctile(cy3_all_molecules(:,12),LowCutoff);%determine high percentile cutoff for dchan spot selection
%dImn=dLowCut;%determine low percentile cutoff for dchan spot selection
%cy3_min_intensity=dI_mean-(dIstd*2.0);%determine low percentile cutoff for dchan spot selection
[dbestMolb, ~,dradXmn,~,~,dradXsd]=curateList(dbestMol,5,3,HighCutoff,LowCutoff);% determine median spot radius X direction
dradXsd=dradXsd*fold_deviation_radius;
cy3_max_x_radius=dradXmn+(dradXsd);
cy3_min_x_radius=dradXmn-(dradXsd/2);
[~, ~,dradYmn,~,~,dradYsd]=curateList(dbestMolb,6,3,HighCutoff,LowCutoff);% determine median spot radius Y direction
dradYsd=dradYsd*fold_deviation_radius;
cy3_max_y_radius=dradYmn+(dradYsd);
cy3_min_y_radius=dradYmn-(dradYsd/2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                Loop Through each segment for spot selection and integration


for FOVind=1:nseg
    


    MovieStart=segments(FOVind,1);%Start point (Frame) to generate average image and ID tehter molecule
    MovieEnd=segments(FOVind,2);%End point (Frame) of trajectory
    segtraj{FOVind}.MovieStart=MovieStart;
    segtraj{FOVind}.MovieEnd=MovieEnd;
    segtraj{FOVind}.SegInfo=segments(segments(:,1)==FOVind,:);%MovieInfo
Cy3ind=1;
Cy5ind=1;
    for n=MovieStart:MovieEnd%sum Cy3 images to ID tethered molecules
        if timeMovie(n,2)==1
            Cy3Time(Cy3ind,1)=timeMovie(n,1);
            Cy3ind=Cy3ind+1;
        elseif timeMovie(n,2)==2
            Cy5Time(Cy5ind,1)=timeMovie(n,1);
            Cy5ind=Cy5ind+1;
        end

    end
    cy5_time=Cy5Time-Cy3Time(1,1);
    cy3_time=Cy3Time-Cy3Time(1,1);

    segtraj{FOVind}.cy5_time=cy5_time;
    segtraj{FOVind}.cy3_time=cy3_time;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%            Plot selection information                   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plotIvsRad_color_select(cy3_all_molecules(:,12),cy3_all_molecules(:,5),cy3_max_intensity,cy3_min_intensity,cy3_max_x_radius,cy3_min_x_radius,'Intensity',...
        'X Radius (Pixels)',25,'b',[0.8500 0.3250 0.0980],'Cy3 Intensity vs Xradius',dataout,outf,FOVind);
    plotIvsRad_color_select(cy3_all_molecules(:,12),cy3_all_molecules(:,6),cy3_max_intensity,cy3_min_intensity,cy3_max_y_radius,cy3_min_y_radius,'Intensity',...
        'Y Radius (Pixels)',25,'b',[0.8500 0.3250 0.0980],'Cy3 Intensity vs Yradius',dataout,outf,FOVind);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%            Save selection information                   %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cy3_thresholds=[cy3_max_intensity cy3_min_intensity; cy3_max_x_radius cy3_min_x_radius; cy3_max_y_radius cy3_min_y_radius];

    segtraj{FOVind}.SpotSelect.cy3_thresholds=cy3_thresholds;%save threshold values
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     Spot selection based on intensity & Radius          %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cy3_curate_spots=segtraj{FOVind}.id_cy3_spots.Curate;
    cy3_curate_spots(:,12)=cy3_curate_spots(:,11)./cy3_scale;%scale intensity of Cy3 molecules
    cy3_curate_spots=cy3_curate_spots(cy3_curate_spots(:,12)>cy3_min_intensity & cy3_curate_spots(:,12)<cy3_max_intensity & cy3_curate_spots(:,5)>=cy3_min_x_radius & cy3_curate_spots(:,5)<=cy3_max_x_radius ...
        & cy3_curate_spots(:,6)>=cy3_min_y_radius & cy3_curate_spots(:,6)<=cy3_max_y_radius,:);
    %keep points with intensity and radius within acceptable range

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Determine the xy coordinates of peaks in cy5 channel             %%
    %%%     based on LSF using local reference spots                      %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    refchannelA=nanogrids.finDpt;
    refchannelB=nanogrids.finApt;
    [cy3_spots_translate_cy5]=mappoints_between_channels(cy3_curate_spots(:,1:2),refchannelA,refchannelB,rectangle,number_Points_to_Fit);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Loop through Cy3 and cy5 spots                                   %%
    %%%    and save pairs based on distance                               %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cy5_id_spots=segtraj{FOVind}.id_cy5_spots.Curate;
    dchDistmat=pdist2(cy5_id_spots(:,1:2),cy3_spots_translate_cy5(:,1:2)); %calculate the pairwise distances between points in Cy3 and Cy5
    [distance_between_cy3_cy5_spots,index_of_spots]=sort(dchDistmat,2);%sort the distances down columns
    c=1;
    cy3_revised_spots=[];
    cy3_translate_cy5_revised_spots=[];
    cy5_revised_spots=[];
    molesLst=[];
    for jj=1:size(cy5_id_spots,1)
        tmp_spot_location=[];
        current_spot_index=find(index_of_spots(:,1)==index_of_spots(jj,1));%Determine if Cy3 spot jj is closest for multiple Cy5 spots
        tmp_spot_location=distance_between_cy3_cy5_spots(current_spot_index,1);%get the distance for each pair
        tmp_spot_location(:,2)=current_spot_index;%get index of Cy5 spot pairs
        [~,index_closest_spot]=min(tmp_spot_location(:,1));%find closest molecule molecular pairs
        tmp_spot_location(:,3)=index_of_spots(index_of_spots(:,1)==index_of_spots(jj,1),1);%Index of molecule in 532nm channel
        tmp_spot_location(1,4)=index_closest_spot;% add Index closest molecule in Cy5 channel
        molesLst{jj}=tmp_spot_location; %save indexing data to array


        if tmp_spot_location(index_closest_spot,1)<max_distance_between_cy3_cy5 &&  tmp_spot_location(index_closest_spot,2)==jj %Keep molecules within a given distance of each other
            cy3_revised_spots(c,1:2)=cy3_curate_spots(index_of_spots(jj,1),1:2);%save Cy3 coordinates that are within cutoff
            cy3_revised_spots(c,3)=index_of_spots(jj,1);%spot number
            cy3_revised_spots(c,4)=tmp_spot_location(1); %distance between pairs
            cy3_revised_spots(c,5)=jj; %spot number in Cy5 channel
            cy3_revised_spots(c,6:7)=cy5_id_spots(jj,1:2);
            cy5_revised_spots(c,1:2)=cy5_id_spots(jj,1:2);
            cy5_revised_spots(c,3)=jj;
            cy5_revised_spots(c,4)=index_of_spots(jj,1);
            cy5_revised_spots(c,5)=tmp_spot_location(index_closest_spot,1);
            c=c+1;
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                                           %%%
    %%%      Determine the Drift correction       %%%
    %%%                                           %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    waitbar(.5,bar_wait,'Calculate Drift');

    mnSgRt=median(cy3_curate_spots(:,11));
    stdSrRT=std(cy3_curate_spots(:,11));
    SigRatio=mnSgRt-stdSrRT;
    drift_info=[];
    [offx,offy,drift_info]=calculate_drift_correction(cy3_curate_spots(:,1:2),reader,MovieStart,MovieEnd,timeMovie,95,border,40,area,sen,SigRatio,0,FOVind,outf);
    close all;
    segtraj{FOVind}.offx = offx;
    segtraj{FOVind}.offy = offy;
    segtraj{FOVind}.drift_info = drift_info;
    maxxdrift = max(offx);
    minxdrift = min(offx);
    maxydrift = max(offy);
    minydrift = min(offy);
  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     Remove points that drift outside FOV during movie             %%
    %%%                                                                   %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Cy5mpidx,~]=find(cy5_revised_spots(:,2)+2*boxwidth+maxydrift<512 & cy5_revised_spots(:,2)-2*boxwidth-minydrift>0 & cy5_revised_spots(:,1)+2*boxwidth...
        +maxxdrift<256 & cy5_revised_spots(:,1)-2*boxwidth-minxdrift>0);
    [Cy3idx,~]=find(cy3_revised_spots(:,2)+2*boxwidth+maxydrift<512 & cy3_revised_spots(:,2)-2*boxwidth-minydrift>0 & cy3_revised_spots(:,1)+ 2*boxwidth+maxxdrift<256 & cy3_revised_spots(:,1)-2*boxwidth-minxdrift>0);

    Pick=intersect(Cy3idx,Cy5mpidx);%find xy coordinates that are consistent between two channels

    cy5_final_spots=round(cy5_revised_spots(Pick(:,1),1:2));
    cy3_final_spots=round(cy3_revised_spots(Pick(:,1),1:2));

    cy3_final_spots(:,3)=1;%Add a boolean for ID of bleach analysis
    cy5_final_spots(:,3)=1;%Add a boolean for ID of bleach analysis

    segtraj{FOVind}.cy3_final_spots=cy3_final_spots; %save coordinates for downstream anaylsis
    segtraj{FOVind}.cy5_final_spots=cy5_final_spots; %save coordinates for downstream anaylsis


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                                           %%%
    %%%   plot rejected vs final spots            %%%
    %%%                                           %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cy3_all_points=segtraj{FOVind}.id_cy3_spots.Centroid;
    cy3_all_points=round(cy3_all_points(:,1:2));
    cy3_saved_pionts=round(cy3_final_spots(:,1:2));
    Lid = ismember(cy3_all_points,cy3_saved_pionts,'rows');
    cy3_rejected_spots=cy3_all_points(Lid==0,:);
    cy3_all_points(Lid==0,9)=0;
    cy3_all_points(Lid==1,9)=1;

    segtraj{FOVind}.cy3_rejected_spots=cy3_rejected_spots;
    [cy5_rejected_spots]=mappoints_between_channels(cy3_rejected_spots(:,1:2),refchannelA,refchannelB,rectangle,number_Points_to_Fit);  %Calculate the postition of badpoints in acceptor channel


    %Plot the position of thether molecules that satisfy intensity value and overlay final thether molecules

    f = figure;
    set(f, 'Visible', 'on'); clf;
    cy3_img_average=segtraj{FOVind}.cy3_img_average; %load average image for comparison
    cy5_img_average=segtraj{FOVind}.cy5_img_average; %load average image for comparison

    set(f,'Position',[477 136 954 827]);
    movegui(gcf,'center');
    ax1=subplot(121);
    imshow(cy3_img_average,'DisplayRange',[mean(quantile(cy3_img_average,0.80)) mean(quantile(cy3_img_average,0.99))]); hold on;
    make_boxplot(cy3_final_spots, cy3_rejected_spots,area,0);
    set(gca,'Ydir','reverse');
    set(gca,'yLim',[0 520]);
    set(gca,'xLim',[0 256]);
    %title('Raw');
    textborder(140, -10, 'Goodpoints', 'green', 'black', 'FontName','Arial Black','FontSize',16,'FontWeight','bold');
    textborder(140, -30, 'BadPoints', 'red', 'black', 'FontName','Arial Black','FontSize',16,'FontWeight','bold');
    textborder(5, 520, ['Cy3'], 'white', 'black', 'FontName','Arial','FontSize',12);
    set(gca,'FontName','Arial Black','FontSize',16,'FontWeight','bold','LineWidth',3)%,'XColor',[0 0 0],'YColor',[0 0 0],'YTick','ZColor',[0 0 0]);
    hold off;
    ax2=subplot(122);
    imshow(cy5_img_average,'DisplayRange',[mean(quantile(cy5_img_average,0.80)) mean(quantile(cy5_img_average,0.99))]); hold on;
    make_boxplot(cy5_final_spots, cy5_rejected_spots,area,0);
    set(gca,'Ydir','reverse');
    set(gca,'yLim',[0 520]);
    set(gca,'xLim',[0 256]);
    %title(' transfomed ');
    textborder(140, -10, 'Goodpoints', 'green', 'black', 'FontName','Arial Black','FontSize',16,'FontWeight','bold');
    textborder(140, -30, 'BadPoints', 'red', 'black', 'FontName','Arial Black','FontSize',16,'FontWeight','bold');
    textborder(5, 520, ['Cy5'], 'white', 'black', 'FontName','Arial','FontSize',12);
    set(gca,'FontName','Arial Black','FontSize',16,'FontWeight','bold','LineWidth',3)%,'XColor',[0 0 0],'YColor',[0 0 0],'YTick','ZColor',[0 0 0]);
    hold off;
    %saveas(gcf, [outf '_Points_overLap_' num2str(segind) '.png']);
    saveas(gcf, fullfile(dataout,[outf num2str(FOVind) '_Points_overLap']), 'png');hold off;
    close all;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Determine nonspecific binding areas or "dark windows" in the tether  %%
    %                       molecule field of view                           %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    waitbar(.75,bar_wait,'ID nonspecific ROI');

    Cy3allB=segtraj{FOVind}.id_cy3_spots.Psf;
    dark_data=[];
    [cy3_drk_roi,dark_data.map,dark_data.fits]=DarkBox_fits(cy3_img_average,border,drkwindow,Cy3allB,boxwidth,cy3_final_spots,FOVind,outf);

   % DarkBox_fits(image,border,minDist,brightspots,boxwidth,FinalPoints,segmentnumber,outname)
    dark_data.points=cy3_drk_roi;
    segtraj{FOVind}.dark_roi_info=dark_data; %save dark info for plotting




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Determine the xy coordinates of darkwindows in Cy5 channel       %%
    %%%     based on LSF using local reference spots                      %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [cy5_drk_roi]=mappoints_between_channels(cy3_drk_roi,refchannelA,refchannelB,rectangle,number_Points_to_Fit);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     Remove points that drift outside FOV during movie             %%
    %%%                                                                   %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   
    [Cy5mpidx,~]=find(cy5_drk_roi(:,2)+2*boxwidth+maxydrift<512 & cy5_drk_roi(:,2)-2*boxwidth-minydrift>0 & cy5_drk_roi(:,1)+2*boxwidth+maxxdrift<256 & cy5_drk_roi(:,1)-2*boxwidth-minxdrift>0);
    [Cy3idx,~]=find(cy3_drk_roi(:,2)+2*boxwidth+maxydrift<512 & cy3_drk_roi(:,2)-2*boxwidth-minydrift>0 & cy3_drk_roi(:,1)+ 2*boxwidth+maxxdrift<256 & cy3_drk_roi(:,1)-2*boxwidth-minxdrift>0);
    Pick=intersect(Cy3idx,Cy5mpidx);%find xy coordinates that are consistent between two channels
    cy5_dark_final_roi=round(cy5_drk_roi(Pick(:,1),1:2));%xy coordinates for dark roi in cy5 channel
    cy3_dark_final_roi=round(cy3_drk_roi(Pick(:,1),1:2));%xy coordinates for dark roi in cy3 channel

    cy3_dark_final_roi(:,3)=1;%Add a boolean for ID of bleach analysis
    cy5_dark_final_roi(:,3)=1;%Add a boolean for ID of bleach analysis
    segtraj{FOVind}.cy5_dark_final_roi=cy5_dark_final_roi; %save coordinates for downstream anaylsis
    segtraj{FOVind}.cy3_dark_final_roi=cy3_dark_final_roi; %save coordinates for downstream anaylsis





    %Plot the donor channel and overlay the xy coordinates for peaks and
    %dark spaces

    f = figure;
    set(f, 'Visible', 'on'); clf;
    set(f,'Position',[477 136 954 827]);
    movegui(gcf,'center');

    subplot(1,2,1);
    imshow(cy3_img_average,'DisplayRange',[mean(quantile(cy3_img_average,0.80)) mean(quantile(cy3_img_average,0.99))]); hold on;
    make_boxplot(cy3_final_spots, cy3_dark_final_roi,area,0);

    set(gca,'Ydir','reverse');
    set(gca,'yLim',[0 520]);
    set(gca,'xLim',[0 256]);
    textborder(140, -10, 'Bright', 'green', 'black', 'FontName','Arial Black','FontSize',16,'FontWeight','bold');
    textborder(140, -30, 'Darkspots', 'red', 'black', 'FontName','Arial Black','FontSize',16,'FontWeight','bold');
    textborder(5, 520, ['Cy3'], 'white', 'black', 'FontName','Arial Black','FontSize',16);
    set(gca,'FontName','Arial Black','FontSize',16,'FontWeight','bold','LineWidth',3)%,'XColor',[0 0 0],'YColor',[0 0 0],'YTick','ZColor',[0 0 0]);
    hold off;

    subplot(1,2,2);

    %Plot the acceptor channel and overlay the xy coordinates for peaks
    imshow(cy5_img_average,'DisplayRange',[mean(quantile(cy5_img_average,0.80)) mean(quantile(cy5_img_average,0.99))]); hold on;
    make_boxplot(cy5_final_spots, cy5_dark_final_roi,area,0);

    textborder(140, -10, 'Brighspots', 'green', 'black', 'FontName','Arial Black','FontSize',16,'FontWeight','bold');
    textborder(140, -30, 'Darkspots', 'red', 'black', 'FontName','Arial Black','FontSize',16,'FontWeight','bold');
    textborder(5, 520, ['Cy5'], 'white', 'black', 'FontName','Arial Black','FontSize',16);
    set(gca,'FontName','Arial Black','FontSize',16,'FontWeight','bold','LineWidth',3)%,'XColor',[0 0 0],'YColor',[0 0 0],'YTick','ZColor',[0 0 0]);
    hold off;
    %saveas(gcf, [outf '_shifted_overlay_' num2str(segind) '.png']);
    saveas(gcf, fullfile(dataout,[outf num2str(FOVind) '_Cy3_Cy5_Bright and Dark Spots']), 'png');hold off;
    close all;

    % Initialize trajectory and background arrays
    Cy5ind=1;
    Cy3ind=1;
    Cy5PTS=[];
    Cy3PTS=[];

    totframe=1;

    for n=MovieStart:MovieEnd%sum Cy3 images to ID tethered molecules
        if timeMovie(n,2)==1
            Cy3PTS(Cy3ind,1)=Cy3ind;
            Cy3ind=Cy3ind+1;
        elseif timeMovie(n,2)==2
            Cy5PTS(Cy5ind,1)=Cy5ind;
            Cy5ind=Cy5ind+1;
        end

        totframe=totframe+1;
    end

    Cy3len=size(Cy3PTS,1);
    Cy5len=size(Cy5PTS,1);
    numpBrght=size(cy3_final_spots,1);
    numpDrk=size(cy3_dark_final_roi,1);


    bxsize=size(sen,1);

    segtrajImg.cy3_images = zeros(bxsize,bxsize,numpBrght,Cy3len);
    segtrajImg.cy3_dark_roi_images= zeros(bxsize,bxsize,numpDrk,Cy3len);
    segtrajImg.cy5_images = zeros(bxsize,bxsize,numpBrght,Cy5len);
    segtrajImg.cy5_dark_roi_images = zeros(bxsize,bxsize,numpDrk,Cy5len);

    segtraj{FOVind}.analyzed=false;

    % include information from psfFit_Image (2D gaussian fit) 
    % output [xpos,ypos,Intensity,Background,sigma_x,sigma_y,angle; quality of fit].
    for jk=1:numpBrght
        bleach_analysis{jk}.changepoint_list=[];
    end
   for jk=1:numpDrk
        bleach_analysis_dark{jk}.changepoint_list=[];
    end
    segtraj{FOVind}.bleach_analysis=bleach_analysis;%empty cell array for bleach analysis
    segtraj{FOVind}.bleach_analysis=bleach_analysis_dark;%empty cell array for bleach analysis
    segtraj{FOVind}.cy3_intensity=zeros(Cy3len,numpBrght);%Intensity of bright spots in cy3 channel
    segtraj{FOVind}.cy3_dark_roi_intensity=zeros(Cy3len,numpDrk);%Intensity of dark spots in cy3 channel
    segtraj{FOVind}.cy3_size=zeros(Cy3len,numpBrght);%size of bright spots in cy3 channel
    segtraj{FOVind}.cy3_dark_roi_size=zeros(Cy3len,numpDrk);%size of dark spots in cy3 channel
    segtraj{FOVind}.cy3_distance=zeros(Cy3len,numpBrght);%distance of bright spots in cy3 channel from tether
    segtraj{FOVind}.cy3_dark_roi_distance=zeros(Cy3len,numpDrk);%size of dark spots in cy3 channel from tether


    segtraj{FOVind}.cy5_intensity=zeros(Cy5len,numpBrght);%Intensity of bright spots in cy5 channel
    segtraj{FOVind}.cy5_dark_roi_intensity=zeros(Cy5len,numpDrk);%Intensity of dark spots in cy5 channel
    segtraj{FOVind}.cy5_size=zeros(Cy5len,numpBrght);%size of bright spots in cy5 channel
    segtraj{FOVind}.cy5_dark_roi_size=zeros(Cy5len,numpDrk);%size of dark spots in cy5 channel
    segtraj{FOVind}.cy5_distance=zeros(Cy5len,numpBrght);%distance of bright spots in cy5 channel from tether
    segtraj{FOVind}.cy5_dark_roi_distance=zeros(Cy5len,numpDrk);%size of dark spots in cy5 channel from tether


    segtraj{FOVind}.cy5_x_position=zeros(Cy5len,numpBrght);%bright spots x-position in cy5 channel
    segtraj{FOVind}.cy5_y_position=zeros(Cy5len,numpBrght);%bright spots y-position in cy5 channel
    segtraj{FOVind}.cy5_dark_roi_x_position=zeros(Cy5len,numpDrk);%dark spots x-position in cy5 channel
    segtraj{FOVind}.cy5_dark_roi_y_position=zeros(Cy5len,numpDrk);%dark spots y-position in cy5 channel







  Cy5ind=1;%this will keep track of Cy5 frame #
    Cy3ind=1;%this will keep track of Cy3 frame #

    
waitbar(1,bar_wait,'Done');
pause(5);
close(bar_wait);

w = waitbar(0,'Integrate Data','Name','fluorescent Molecule integration');


    for j= MovieStart:MovieEnd%loop through movie

        currframe = bfGetPlane(reader,j);
        dch = double(currframe(:,1:256)); % donor half of the image
        ach = double(currframe(:,257:end)); % acceptor half of the image

        waitbar((j-MovieStart+1)/totframe, w, ['Integrating spot intensities in segment ' ...
            num2str(FOVind) ': Frame ' num2str(j-MovieStart+1)...
            ' of ' num2str(totframe)]);

        %Determine the position of each peak based on the drift correction

        current_cy3_position= [cy3_final_spots(:,1)-offx(j,1) cy3_final_spots(:,2) - offy(j,1)];
        current_cy5_position=[cy5_final_spots(:,1) - offx(j,1) cy5_final_spots(:,2) - offy(j,1)];
        %Determine the position of each dark spot based on the drift correction
        current_darkroi_cy3_position= [cy3_dark_final_roi(:,1) - offx(j,1) cy3_dark_final_roi(:,2) - offy(j,1)];
        current_darkroi_cy5_position=[cy5_dark_final_roi(:,1) - offx(j,1) cy5_dark_final_roi(:,2) - offy(j,1)];

        if timeMovie(j,2)==1 % if this is an Cy3ex frame

            %Integrate spot intensity for molecules in cy3 channel
            %Model spots as 2D gaussian for the shape and position of the spots
            [cy3_fit_values,cy3_molecule_image]=integrate_cy3_molecules(dch,current_cy3_position,sen);
            [cy3_dark_roi_fit_values,cy3_dark_roi_molecule_image]=integrate_cy3_molecules(dch,current_darkroi_cy3_position,sen);


            segtraj{FOVind}.cy3_intensity(Cy3ind,:)=(cy3_fit_values(:,3))';%Intensity of bright spots in cy3 channel
            segtraj{FOVind}.cy3_dark_roi_intensity(Cy3ind,:)=(cy3_dark_roi_fit_values(:,3))';%Intensity of dark spots in cy3 channel
            segtraj{FOVind}.cy3_size(Cy3ind,:)=(cy3_fit_values(:,4))';%size of bright spots in cy3 channel
            segtraj{FOVind}.cy3_dark_roi_size(Cy3ind,:)=(cy3_dark_roi_fit_values(:,4))';%size of dark spots in cy3 channel
            segtraj{FOVind}.cy3_distance(Cy3ind,:)=(cy3_fit_values(:,5))';%distance of bright spots in cy3 channel from tether
            segtraj{FOVind}.cy3_dark_roi_distance(Cy3ind,:)=(cy3_dark_roi_fit_values(:,5))';%size of dark spots in cy3 channel from tether

            segtrajImg.cy3_images(:,:,:,Cy3ind)=cy3_molecule_image;
            segtrajImg.cy3_dark_roi_images(:,:,:,Cy3ind)=cy3_dark_roi_molecule_image;

            Cy3ind=Cy3ind+1;

        elseif timeMovie(j,2)==2 % if this is a Cy5ex frame
            %Integrate spot intensity for molecules in cy3 channel
            %Model spots as 2D gaussian for the shape and position of the spots
            [cy5_fit_values,cy5_image]=integrate_cy5_molecules(ach,current_cy5_position,border,pixel_limit_from_cy3,BGcut,BoxDefineLocalBkgrd,area,sen);
            [cy5_dark_roi_fit_values,cy5_dark_roi_molecule_image]=integrate_cy5_molecules(ach,current_darkroi_cy5_position,border,pixel_limit_from_cy3,BGcut,BoxDefineLocalBkgrd,area,sen);

            segtraj{FOVind}.cy5_intensity(Cy5ind,:)=(cy5_fit_values(:,3))';%Intensity of bright spots in cy5 channel
            segtraj{FOVind}.cy5_dark_roi_intensity(Cy5ind,:)=(cy5_dark_roi_fit_values(:,3))';%Intensity of dark spots in cy5 channel
            segtraj{FOVind}.cy5_size(Cy5ind,:)=(cy5_fit_values(:,4))';%size of bright spots in cy5 channel
            segtraj{FOVind}.cy5_dark_roi_size(Cy5ind,:)=(cy5_dark_roi_fit_values(:,4))';%size of dark spots in cy5 channel
            segtraj{FOVind}.cy5_distance(Cy5ind,:)=(cy5_fit_values(:,5))';%distance of bright spots in cy5 channel from tether
            segtraj{FOVind}.cy5_dark_roi_distance(Cy5ind,:)=(cy5_dark_roi_fit_values(:,5))';%size of dark spots in cy5 channel from tether

            segtraj{FOVind}.cy5_x_position(Cy5ind,:)=(cy5_fit_values(:,1))';%bright spots x-position in cy5 channel
            segtraj{FOVind}.cy5_y_position(Cy5ind,:)=(cy5_fit_values(:,2))';%bright spots y-position in cy5 channel
            segtraj{FOVind}.cy5_dark_roi_x_position(Cy5ind,:)=(cy5_dark_roi_fit_values(:,1))';%dark spots x-position in cy5 channel
            segtraj{FOVind}.cy5_dark_roi_y_position(Cy5ind,:)=(cy5_dark_roi_fit_values(:,2))';%dark spots y-position in cy5 channel

            segtrajImg.cy5_images(:,:,:,Cy5ind)=cy5_image;
            segtrajImg.cy5_dark_roi_images(:,:,:,Cy5ind)=cy5_dark_roi_molecule_image;

            Cy5ind=Cy5ind+1;

        end

    end

    close all;

end
save(fullfile(dataout,[outf '_Mol_image.mat']),'segtrajImg','-v7.3','-nocompression' );
st = segtraj; 
close(w)
Done = waitbar(0,'Saving output');
save(fullfile(dataout,[outf '.mat']),'st','-v7.3','-nocompression' );
pause(30);
close(Done);

end


