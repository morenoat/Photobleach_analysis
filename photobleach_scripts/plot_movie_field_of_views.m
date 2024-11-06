function [start_stop_FOV ]= plot_movie_field_of_views(totTime,Cy3_1stFrame,fcorrmatLeftSide,fcorrmatRightSide,dataout)
%%% This function will plot the intensity correlation between periods of Cy3 emission
%%% In the figure the dashed black horizontal line denotes a "threshold" for change in Cy3 intensity between periods, the 
%%% vertical green line indicates first cy3 frame of movie and vertical red
%%% line indicates last frame of movie
%%% A dialog box is used to define the number of FOV, initial Cy3 frame, and last frame of movie
%%% After assigning FOV number, position the cursor at the start of first
%%% FOV and right click on the mouse. Next position the cursor at the end of the FOV and right click with the cursor
%%%A new figure will appear indicating start of FOV with green vertical
%%%lines and stop of FOV as red vertical lines if the assignment is correct
%%%click true, if it is wrong click false and the above process will repeat
corrthresh=0.8;% calculate threshold for plotting horizontal line as a guide for abrupt change in
corrplot = quantile(fcorrmatLeftSide,corrthresh);
corr_cy5_plot = quantile(fcorrmatRightSide,corrthresh);
% correlation indicating a new field of view
Cy3time=totTime(totTime(:,2)==1,1);
Cy5time=totTime(totTime(:,2)==2,1);

marksGood=0;
while marksGood==0

    f = figure;
    set(f, 'Visible', 'on');
    clf; hold off; plot(Cy3time(1:size(fcorrmatLeftSide,2))',fcorrmatLeftSide,'b-o','LineWidth',3);
    hold on; plot(Cy3time(1:size(fcorrmatLeftSide,2))',ones(1,numel(fcorrmatLeftSide))*corrplot,'k--','LineWidth',3);
    plot([totTime(Cy3_1stFrame,1) totTime(Cy3_1stFrame,1)],[0 1],'g','Linewidth',3);
    plot([totTime(end,1)-2 totTime(end,1)-2],[0 1],'r','Linewidth',3);
    set(gcf,'Position', [575 141 793 641]);
    movegui('center');
    set(gca, 'ylim', [0.0 1.5]);
    xlabel('Time sec');
    ylabel('Frame correlation');
    title('Right click to define each FOV start and stop');
    set(gca,'FontName','Arial Black','FontSize',16,'FontWeight','bold','LineWidth',3);
    box off;
    prompt = {'How many FOVs are there?','start time first Segment (s)','End time last segment (s)'};
    dlg_title = 'Define number of FOV';
    dims = [1 35];
    defaultans = {'1',num2str(totTime(Cy3_1stFrame,1)),num2str(totTime(end,1))};
    opts.Resize = 'on';
    answer = inputdlg(prompt,dlg_title,dims,defaultans,opts);


    [xSeg,~] = ginput(str2double(answer{1})*2);
    id_fov=[];
    for hh =1:size(xSeg,1)
        [~,Cy3timeMin]=min(abs(xSeg(hh,1)-Cy3time(:,1)));
        [~,Cy5timeMin]=min(abs(xSeg(hh,1)-Cy5time(:,1)));
        [~,timeMin]=min(abs(Cy3time(Cy3timeMin,1)-totTime(:,1)));
        id_fov(hh,1)=totTime(timeMin,1);
        id_fov(hh,2)=timeMin;
        id_fov(hh,3)=Cy3time(Cy3timeMin,1);
        id_fov(hh,4)=Cy3timeMin;
        id_fov(hh,5)=Cy5time(Cy5timeMin,1);
        id_fov(hh,6)=Cy5timeMin;
    end

    f = figure;
    set(f, 'Visible', 'on'); clf; hold off; plot(Cy3time(1:size(fcorrmatLeftSide,2))',fcorrmatLeftSide,'b-o','LineWidth',3);
    hold on; plot(Cy3time(1:size(fcorrmatLeftSide,2))',ones(1,numel(fcorrmatLeftSide))*corrplot,'k--','LineWidth',3);
    plot([totTime(Cy3_1stFrame,1) totTime(Cy3_1stFrame,1)],[0 1],'r','Linewidth',3);
    plot([totTime(end,1)-2 totTime(end,1)-2],[0 1],'r','Linewidth',3);
    for jk=1:2:size(id_fov,1)-1
        plot([id_fov(jk,1) id_fov(jk,1)],[0 1],'g','Linewidth',3);
        plot([id_fov(jk+1,1) id_fov(jk+1,1)],[0 1],'r','Linewidth',3);
    end
    set(gcf,'Position', [575 141 793 641]);
    movegui('center');
    set(gca, 'ylim', [0.0 1.5]);
    xlabel('Time sec');
    ylabel('Frame correlation');
    title('Correlation between Cy3 frames');
    set(gca,'FontName','Arial Black','FontSize',16,'FontWeight','bold','LineWidth',3);
    box off;
    % Create a dialog box with checkbox selection
    choice = questdlg('Are the selections correct:', 'Repeat FOV Selections?', 'True', 'False','True');

    % Process the selected option
    switch choice
        case 'True'
            disp('True selected');
            marksGood=1;
            plot(-.1:1:1)
            saveas(gcf, fullfile(dataout,'Correlation_Cy3_channel'), 'png');
            close;
        case 'False'
            disp('False selected');
            marksGood=0;
    end



end

saveas(gcf, fullfile(dataout,'Correlation_Cy3_channel'), 'png');

f = figure;
set(f, 'Visible', 'on'); clf; hold off; plot(Cy3time(1:size(fcorrmatRightSide,2))',fcorrmatRightSide,'b-','LineWidth',3);
set(gcf,'Position', [575 141 793 641]);
movegui('center');
set(gca, 'ylim', [0.2 1.2]);
% overlay correlation threshold for FOV segmentation:
hold on; plot(Cy3time(1:size(fcorrmatRightSide,2))',ones(1,numel(fcorrmatRightSide))*corr_cy5_plot,'k--','LineWidth',3);
plot([totTime(end,1)-2 totTime(end,1)-2],[0 1],'r','Linewidth',3);
for jk=1:2:size(id_fov,1)-1
    plot([id_fov(jk,1) id_fov(jk,1)],[0 1],'g','Linewidth',3);
    plot([id_fov(jk+1,1) id_fov(jk+1,1)],[0 1],'r','Linewidth',3);
end
set(gcf,'Position', [575 141 793 641]);
movegui('center');
set(gca, 'ylim', [0.0 1.5]);
xlabel('Time sec');
ylabel('Frame correlation');
title('Correlation between Cy5 frames');
set(gca,'FontName','Arial Black','FontSize',16,'FontWeight','bold','LineWidth',3);
box off;

saveas(gcf, fullfile(dataout,'Correlation_Cy5_channel'), 'png');
movieTime=totTime;






start_stop_FOV=[];
ck=1;
for jk=1:2:size(id_fov,1)-1
    start_stop_FOV(ck,1:2)=id_fov(jk:jk+1,2);%index of FOV start/stop totTime
    start_stop_FOV(ck,3:4)=id_fov(jk:jk+1,1);%time of FOV start/stop totTime
    start_stop_FOV(ck,5:6)=id_fov(jk:jk+1,4);%index of FOV start/stop Cy3Time
    start_stop_FOV(ck,7:8)=id_fov(jk:jk+1,3);%time of FOV start/stop Cy3Time
    start_stop_FOV(ck,9:10)=id_fov(jk:jk+1,6);%index of FOV start/stop Cy5Time
    start_stop_FOV(ck,11:12)=id_fov(jk:jk+1,5);%time of FOV start/stop Cy5Time
    ck=ck+1;
end

close all;
save(fullfile(dataout,'Segments.mat'),'start_stop_FOV');
save(fullfile(dataout,'Movie_time'),'movieTime');
end
