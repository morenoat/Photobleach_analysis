function [start_stop_FOV] = define_movie_segments(segment_output_name,Name_integrated_output,NameDataFile,reference,varargin)
%define_movie_segments(parms_output_folder_name,output_name_data_integration,input_name_data_integration,name_reference_file);

% this is written as a function so that it analyzes the file that you want
% this way, you can analyze different files on different cores
% parallelization


dataout=[pwd '/' segment_output_name];
loc_file=fullfile(dataout,'parms_in.txt');
filename=NameDataFile;
start_stop_FOV=[];
inargs = importdata(loc_file);  % This will give trouble if your file is strictly numeric (it will import it as an array)
Cy3_1stFrame= inargs.data(ismember(inargs.textdata, 'starting_frame')); % Frame number to start analyzing movie 
last_frame= inargs.data(ismember(inargs.textdata, 'last_movie_frame'));% Frame number to end analyzing movie 
number_cy3_excitation_frames= inargs.data(ismember(inargs.textdata, 'number_cy3_excitation_frames'));%number of Cy3 excitation frames
number_cy5_excitation_frames= inargs.data(ismember(inargs.textdata, 'number_cy5_excitation_frames'));%number of Cy5 excitation frames
exposuret=inargs.data(ismember(inargs.textdata, 'exposuret')); %Time of Exposure
delay=inargs.data(ismember(inargs.textdata, 'delay'));% delay between Exposure
border = inargs.data(ismember(inargs.textdata, 'border')); % Width of murky boundary between channels to exclude

if (~isempty(varargin))
    argnum=1;
    while argnum <= numel(varargin)
        switch varargin{argnum}
            case {'Segments'}
                fldr=varargin{argnum+1};
                SegInput=[pwd '/' fldr];
                loc_Seg=fullfile(SegInput,'Segments.mat');
                load(loc_Seg,'start_stop_FOV') ;  % This list the segments
                loc_Movie=fullfile(SegInput,'Movie_time.mat');
                load(loc_Movie,'movieTime');  % This list the time of movie
                argnum=argnum+2;
            otherwise
                error(['Invalid optional argument, ', ...
                    varargin{c}]);
        end
    end
end

if (isempty(start_stop_FOV))
    
    period=number_cy3_excitation_frames+number_cy5_excitation_frames;%defines the period for alternating excitation
    reader = bfGetReader(filename);%read in movie for analysis
    frameID=[];%create an array of alternating excitation period
    frameID(1:number_cy3_excitation_frames,1)=1;% set Cy3 emission in the array as 1
    frameID(number_cy3_excitation_frames+1:number_cy3_excitation_frames+number_cy5_excitation_frames,1)=2;% set Cy3 emission in the array as 2
    nframes = reader.getImageCount;%Determine number of frames in movie
    loopnum=ceil(nframes/period);%determine the number of periods in the movie. 
    cycleframes=repmat(frameID,loopnum,1);%create an array ID every frame as Cy3 or Cy5 emission
    [fcorrmatLeftSide,fcorrmatRightSide] = frame_correlation(reader,period,Cy3_1stFrame,border);%calculate the intensity correlation between periods use this to ID FOVs 
    totTime=(1:1:last_frame).*(exposuret+delay);% create an array of movie time in seconds
    totTime=totTime';
    totTime(Cy3_1stFrame:end,2)=cycleframes(1:size(totTime,1)-Cy3_1stFrame+1,1);%label time array with Cy3 or Cy5 emission
    fcorrmatLeftSide = (fcorrmatLeftSide-min(fcorrmatLeftSide))./(max(fcorrmatLeftSide) - min(fcorrmatLeftSide));%normalize period correlation for Cy3 channel
    fcorrmatRightSide = (fcorrmatRightSide-min(fcorrmatRightSide))./(max(fcorrmatRightSide) - min(fcorrmatRightSide));%normalize period correlation for Cy5 channel

    [start_stop_FOV ]= plot_movie_field_of_views(totTime,Cy3_1stFrame,fcorrmatLeftSide,fcorrmatRightSide,dataout);% plot correlation vs time and assign start and stop of FOVs
 


end


warning('off','all');
warning
segnum=start_stop_FOV;
timeMov=movieTime;
spot_selection_integration_frame_by_frame(filename,loc_file,Name_integrated_output ,reference,'segments', segnum,'MovieTime',timeMov);


end
