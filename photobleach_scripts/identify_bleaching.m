function [data]=identify_bleaching(data,parameters)

scale_value= parameters.scale_value;%set max intensity to ID final photobleaching step
smoothing=parameters.smoothing;%how much to smooth trajectories to determine max intensity this will get rid of outlier points
minimum_intensity=parameters.minimum_intensity;%0.2 After inital background subtration set limit for if molecules have been corrected
fractional_change=parameters.fractional_change; % 0.3 minimum for fractional change in cy5 signal to be considered photobleaching event
relative_to_max_intensity=parameters.relative_to_max_intensity;%2threhold to ID state changes relative to max intensity
move_median_sample=parameters.move_median_sample;%Number of points to calculate moving mean
minimum_points=parameters.minimum_points;%remove transient states
replace_negative= parameters.replace_negative ;%value to replace negative intensity values
dark_roi=parameters.dark_roi;%analyze the dark roi (nonspecific bound) for photobleaching
cy5_time=data.cy5_time;
mol_list=data.cy5_final_spots(:,3);
if dark_roi==false
    trajectories=data.cy5_intensity;
    cy3_data=data.cy3_intensity;
else
    trajectories=data.cy5_dark_roi_intensity;
    cy3_data=data.cy3_dark_roi_intensity;
end

mnCy5=movmedian(trajectories,move_median_sample,1);%smooth data for photobleaching analysis
mnCy5b=movmedian(trajectories,smoothing,1);%more aggressive smothing to determine max intensity
scale_cy5_intensity=[];
mean_cy5_max=max(mnCy5b,[],1);
for jh=1:size(mnCy5,2)
    scale_cy5_intensity(:,jh)=mnCy5(:,jh)./(mean_cy5_max(1,jh)/scale_value);%scale the data based on Intensity_limit
end


mnCy3=movmedian(cy3_data,move_median_sample,1);%smooth data for photobleaching analysis
mnCy3b=movmedian(cy3_data,smoothing,1);%more aggressive smothing to determine max intensity
smooth_cy3_intensity=[];
scale_cy3_intensity=[];
mean_cy3_max=max(mnCy3b,[],1);
for jh=1:size(mnCy3,2)
    smooth_cy3_intensity(:,jh)=mnCy3(:,jh)./(mean_cy3_max(1,jh)/scale_value);%scale the data based on Intensity_limit
    scale_cy3_intensity(:,jh)=cy3_data(:,jh)./(mean_cy3_max(1,jh)/scale_value);%scale the data based on Intensity_limit
end



intensity_data = flip(scale_cy5_intensity,1);%Flip Trajectory to ID last photobleaching event
num_trajectories = size(intensity_data,2);%number of trajectories
% Find the point of loss of signal intensity for each trajectory
IDCHANGEA=[];%keep alist of changepoints and how much the intensity has changed
intensity_data_b=[];

for i = 1:num_trajectories
    point_of_change=[];
    point_of_change_b=[];
    % Replace this with your actual logic to find the loss point
    [binary_of_change,smoothed_intensity,~]=ischange(intensity_data(:, i),'mean','Threshold',scale_value);%
    %ID changepoints based on mean intensity change

    if any(binary_of_change>0)
        point_of_change=find(binary_of_change==1);%ID pos of Change
        point_of_change=[1;point_of_change];%list of starting points for intensity states
        point_of_change(1:end-1,2)=point_of_change(2:end,1);%previous changepoint
        point_of_change(end,2)=size(smoothed_intensity,1);%add end of trajectory
        point_of_change(end,:)=[];%now delete the end
        point_of_change(:,3)=smoothed_intensity(point_of_change(:,2)-1,1);%Intensity of previous state
        point_of_change(:,4)=smoothed_intensity(point_of_change(:,2),1);%Intensity of new state
        point_of_change(:,5)=(point_of_change(:,4)-point_of_change(:,3));%difference in intensity between states
        point_of_change(point_of_change(:,5)<0,:)=[];%remove changepoints that represent intensity increase
        if ~isempty(point_of_change)
            point_of_change(:,5)=(point_of_change(:,4)-point_of_change(1,3));%difference in intensity relative to background
            point_of_change(:,6)=(point_of_change(:,5)./point_of_change(1,5));%normalize subtracted intensity to pre-background intensity
            IDCHANGEA{i}.point_of_change_initial=point_of_change;
            intensity_data_b(:,i)=intensity_data(:, i)-point_of_change(1,3);%difference in intensity relative to background
            intensity_data_b(:,i)=intensity_data_b(:,i)./point_of_change(1,5);%normalize subtracted intensity to pre-background intensity
            smoothed_intensity=smoothed_intensity(:, 1)-point_of_change(1,3);%difference in intensity relative to background
            smoothed_intensity(:,1)=smoothed_intensity(:,1)./point_of_change(1,5);%normalize subtracted intensity to pre-background intensity
            for bf=1:size(point_of_change,1)
                point_of_change(bf,7)=mean(intensity_data_b(point_of_change(bf,1):point_of_change(bf,2),i));%mean intensity of intial state
            end
            point_of_change(1:end-1,8)=point_of_change(2:end,7);%mean intensity of next state
            point_of_change(end,8)=mean(intensity_data_b(point_of_change(end,2):end,i));%mean intensity of first state

            %%% Now reverse the order of the changepoints

            point_of_change_b=[1;(size(smoothed_intensity,1)-point_of_change(end:-1:1,2)+1)];
            point_of_change_b(1:end-1,2)=point_of_change_b(2:end,1);%previous changepoint
            point_of_change_b(end,2)=size(smoothed_intensity,1);%add end of trajectory
            point_of_change_b(end,:)=[];%now delete the end
            point_of_change_b(:,3)=point_of_change(end:-1:1,8);
            point_of_change_b(:,4)=point_of_change(end:-1:1,7);
            point_of_change_b(:,8)=1;%Denote bleaching events
        else
            intensity_data_b(:,i)=intensity_data(:, i);
            point_of_change_b(1,1)=1;
            point_of_change_b(1,2)=size(smoothed_intensity,1);%add end of trajectory
            point_of_change_b(1,3)=mean(intensity_data_b(1:end,i));%mean intensity of first state
            point_of_change_b(1,4)=mean(intensity_data_b(1:end,i));%mean intensity of first state
            point_of_change_b(1,8)=0;%mean intensity of first state
        end
    else
        intensity_data_b(:,i)=intensity_data(:, i);
        point_of_change_b(1,1)=1;
        point_of_change_b(1,2)=size(smoothed_intensity,1);%add end of trajectory
        point_of_change_b(1,3)=mean(intensity_data_b(1:end,i));%mean intensity of first state
        point_of_change_b(1,4)=mean(intensity_data_b(1:end,i));%mean intensity of first state
        point_of_change_b(1,8)=0;%mean intensity of first state

    end
    IDCHANGEA{i}.changepoint_list=point_of_change_b;
    IDCHANGEA{i}.point_of_change_initial=point_of_change;
    IDCHANGEA{i}.smoothed_intensity=flip(smoothed_intensity,1);
    IDCHANGEA{i}.normalized_intensity=flip(intensity_data_b(:,i),1);

end
scale_cy5_intensity= flip(intensity_data_b,1);
cy5_fit_intensity=[];

IDCHANGE=[];
for hh=1:size(scale_cy5_intensity,2) %Normalize trajectory intensity to Itensity preceding last bleaching event and ID changepts. only keep ChangePts preceding bleach event
    point_of_change=[];
    more_postive_changes=1;%boolean for removing changepoints with increasing intensity
    more_transient_changes=1;%boolean for removing transient changepoints
    normalized_intensity=scale_cy5_intensity(:,hh);
    max_intensity=max(normalized_intensity);
    [binary_of_change,smoothed_intensity,S2]=ischange(normalized_intensity,'mean','Threshold',max_intensity*relative_to_max_intensity);%ID change points and Signal of States and Variance
    IDCHANGE{hh}.moving_stdev=S2;
    cy5_fit_intensity(:,hh)=smoothed_intensity;
    if any(binary_of_change>0)
        point_of_change=find(binary_of_change==1);%ID pos of Change
        point_of_change=[1;point_of_change];
        point_of_change(1:end-1,2)=point_of_change(2:end,1);%previous changepoint
        point_of_change(end,2)=size(normalized_intensity,1);%add end of trajectory
        point_of_change(end,:)=[];%now delete the end
        point_of_change(:,3)=smoothed_intensity(point_of_change(:,2)-1,1);%Intensity of previous state
        point_of_change(:,4)=smoothed_intensity(point_of_change(:,2),1);%Intensity of new state
        %         point_of_change(point_of_change(:,4)>point_of_change(:,3),:)=[];%remove changepoints that increase
        IDCHANGE{hh}.initial_points=point_of_change;
        while more_transient_changes==1 %this loop will remove transient change points
            if ~isempty(point_of_change)
                point_of_change(:,1)=[1;point_of_change(1:end-1,2)];
                point_of_change(1:end-1,2)=point_of_change(2:end,1);%previous changepoint
                transient_points=[];
                [transient_points,~]=find(point_of_change(:,2)-point_of_change(:,1)<minimum_points);
                if ~isempty(transient_points)
                    point_of_change(transient_points(1,1),:)=[];
                else
                    more_transient_changes=0; %exit out of loop
                end

            else
                more_transient_changes=0; %exit out of loop
            end
        end
        while more_postive_changes==1
            if ~isempty(point_of_change)
                point_of_change(:,1)=[1;point_of_change(1:end-1,2)];
                for ww =1:size(point_of_change,1)
                    point_of_change(ww,3)=mean(normalized_intensity(point_of_change(ww,1):point_of_change(ww,2),1));%Intensity of previous state
                end
                point_of_change(1:end-1,4)=point_of_change(2:end,3);%Intensity of the next state
                point_of_change(end,4)=mean(normalized_intensity(point_of_change(end,2):end,1));%Intensity of last state
                point_of_change(point_of_change(:,3)<0,3)=replace_negative;%replace negative values with min value
                point_of_change(point_of_change(:,4)<0,4)=replace_negative;%replace negative values with min value
                point_of_change(:,5)=1-(point_of_change(:,4)./point_of_change(:,3));%calculate fractional change between states
                transition_points_to_remove=[];
                [transition_points_to_remove,~]=find(point_of_change(:,5)<fractional_change...
                    | point_of_change(:,3)<minimum_intensity);
                if ~isempty(transition_points_to_remove)
                    point_of_change(transition_points_to_remove(1,1),:)=[];%remove changepoints that increase
                else
                    more_postive_changes=0; %exit from loop
                end

            else
                more_postive_changes=0; %exit from loop

            end
        end

        %%%%NEW Section
        if size(point_of_change,1)>0
            [transition_points_to_remove,~]=find( point_of_change(:,4)<minimum_intensity,1);
            if transition_points_to_remove<size(point_of_change,1)
                point_of_change(transition_points_to_remove+1:end,:)=[];%remove transistions after initial photobleach
            end


        end
        if size(point_of_change,1)<1
            point_of_change(1,1)=1;
            point_of_change(1,2)=size(normalized_intensity,1);
            point_of_change(1,3)=mean(normalized_intensity(point_of_change(1,1):point_of_change(1,2),1));%median State Intensity
            point_of_change(1,4)=mean(normalized_intensity(point_of_change(1,1):end,1));%median post State Intensity
            point_of_change(1,8)=0;%Denote no bleaching events
        else
            point_of_change(:,8)=1;%Denote bleaching events
        end
    else
        point_of_change(1,1)=1;
        point_of_change(1,2)=size(normalized_intensity,1);
        point_of_change(1,3)=mean(normalized_intensity(point_of_change(1,1):point_of_change(1,2),1));%median State Intensity
        point_of_change(1,4)=mean(normalized_intensity(point_of_change(1,1):end,1));%median post State Intensity
        point_of_change(1,8)=0;%Denote no bleaching events
    end

    if mol_list(hh,1)==0
        point_of_change=[];
         point_of_change(1,1)=1;
        point_of_change(1,2)=size(normalized_intensity,1);
        point_of_change(1,3)=mean(normalized_intensity(point_of_change(1,1):point_of_change(1,2),1));%median State Intensity
        point_of_change(1,4)=mean(normalized_intensity(point_of_change(1,1):end,1));%median post State Intensity
        point_of_change(1,8)=0;%Denote no bleaching events
    end


    point_of_change(:,5)=cy5_time(point_of_change(:,1),1);%convert to time
    point_of_change(:,6)=cy5_time(point_of_change(:,2),1);%convert to time
    point_of_change(:,7)=hh;

   % bleached_intensity=vertcat(bleached_intensity,point_of_change(end,4));%make a list of background intensity
   % last_molecule_intensity=vertcat(last_molecule_intensity,point_of_change(end,3));%make a list of background intensity
   % initial_molecule_intensity=vertcat(initial_molecule_intensity,point_of_change(1,3));%make a list of initial intensity
    IDCHANGE{hh}.last_bleach=point_of_change(end,:);
    IDCHANGE{hh}.changepoint_list=point_of_change;



end




change_point_list=[];%Number of states per Trajectory
c=1;
all_molecules=[];
molecules_that_bleach=[];
molecules_that_bleach_cy5_intensity=[];
molecules_that_bleach_cy3_intensity=[];
number_of_molecules=size(scale_cy5_intensity,2);
vv=1;
for tx=1:size(IDCHANGE,2)
    tmpLst=[];
    tmpLst=IDCHANGE{tx}.changepoint_list(IDCHANGE{tx}.changepoint_list(:,8)==1,:);
    mrkChg=tmpLst;
    if size(mrkChg,1)>0
        mrkChg(:,7)=tx;
        mrkChg(:,8)=1:1:size(mrkChg,1);
        change_point_list(c,1)=size(mrkChg,1);
        change_point_list(c,3)=mrkChg(end,3);%Intensity of second to last state
        change_point_list(c,4)=mrkChg(end,4);%Intensity of the bleach state
        change_point_list(c,2)=mrkChg(1,3);%Intensity of the initial state
        all_molecules=vertcat(all_molecules,mrkChg);
        molecules_that_bleach{vv}=IDCHANGE{tx};
        molecules_that_bleach_cy5_intensity(:,vv)=scale_cy5_intensity(:,tx);
        molecules_that_bleach_cy3_intensity(:,vv)=smooth_cy3_intensity(:,tx);
        vv=vv+1;
    else
        change_point_list(c,1)=0;
        change_point_list(c,2)=scale_cy5_intensity(1,tx);
        change_point_list(c,3)=scale_cy5_intensity(1,tx);
        change_point_list(c,4)=scale_cy5_intensity(1,tx);
    end
    c=c+1;
end

all_molecules(:,6)=all_molecules(:,2)-all_molecules(:,3);
list_bleaching_events=change_point_list(change_point_list(:,1)>0,:);
tabulate_bleaching=[];
tabulate_bleaching.all_molecules=all_molecules;
tabulate_bleaching.list_bleaching_events=list_bleaching_events;
if dark_roi==false
    data.cy5_tabulate_bleaching=tabulate_bleaching;
    data.cy5_intensity_scale=scale_cy5_intensity;
    data.cy5_fit_intensity=cy5_fit_intensity;
    data.cy3_intensity_smooth=smooth_cy3_intensity;
    data.cy3_intensity_scale=scale_cy3_intensity;
    data.bleach_analysis=IDCHANGE;
else
    data.cy5_dark_roi_tabulate_bleaching=tabulate_bleaching;
    data.cy5_dark_roi_intensity_scale=scale_cy5_intensity;
    data.cy5_dark_roi_fit_intensity=cy5_fit_intensity;
    data.cy3_dark_roi_intensity_smooth=smooth_cy3_intensity;
    data.cy3_dark_roi_intensity_scale=scale_cy3_intensity;
    data.bleach_analysis_dark=IDCHANGE;
end



end
