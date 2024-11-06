function plotIvsRad_color_select(Intensity,Radius,max_intensity,min_intensity,max_radius,min_radius,Xaxis,Yaxis,spot_size,color_all_points,color_select_points,plotName,dataout,outf,segment_number)
    
    all_points=[(Intensity),Radius];
    
    select_points=all_points(all_points(:,1)>=min_intensity & all_points(:,1)<=max_intensity & all_points(:,2)>=min_radius & all_points(:,2)<=max_radius,:);
    f = figure;
    set(f, 'Visible', 'on'); clf;
    set(f,'Position',[666 188 950 581])
    scatter(all_points(:,1),all_points(:,2),spot_size,color_all_points); hold on;
    scatter(select_points(:,1),select_points(:,2),spot_size,color_select_points); hold on;
    X = prctile(Intensity,99.0);
    Y = prctile(Radius,70);
    %line([max_intensity max_intensity],[0, 10],'LineWidth',3,'Color','red','LineStyle','--');%plot line of 98% 
    %line([min_intensity  min_intensity],[0, 10],'LineWidth',3,'Color','red','LineStyle','--');%plot line of 98% 
    % hline = refline([0  max_radius]);
    % hline.Color = 'r';
    % hline.LineStyle = '--';
    % hline.LineWidth=3;
    % lwline = refline([0  min_radius]);%'Color','red','LineStyle','--');%plot line of 98% 
    % lwline.Color = 'r';
    % lwline.LineStyle = '--';
    lwline.LineWidth=3;
    ylabel(Yaxis)
    xlabel(Xaxis)
   
    set(gca,'Ylim',[0 6])
    %set(gca,'Xlim',[-.01 X])
    title(plotName);
    box off;
    set(gca,'FontName','Arial Black','FontSize',16,'FontWeight','bold','LineWidth',3);
    saveas(gcf, fullfile(dataout,[outf '_' num2str(segment_number) '_' plotName]), 'png');hold off;
    close all;
    
end

