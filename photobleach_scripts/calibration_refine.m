function []=calibration_refine(inputA,nangridsImageNumber,rectangle,number_Points_to_Fit,CutoffDistance,plot_each_image,outf,path_to_data,varargin)
nanogrids=importdata(inputA);
%rectangle = size of box to get reference data points
numngs=nangridsImageNumber;%number of reference images
minNumPt=number_Points_to_Fit;%Number of reference necessary for least squares fitting
minDistanceBetweenPointsB=CutoffDistance;%maximum distance between translated and FOV point
minDistanceBetweenPointsA=CutoffDistance;%maximum distance between translated and FOV point
finalRef=[];
plot_final=false;
plot_nanogrid_final=false;
if (~isempty(varargin))
    argnum=1;
    while argnum <= numel(varargin)
        switch varargin{argnum}
            case {'Last_cycle'} %do you want to plot reference set overlap
                plot_final=varargin{argnum+1};
                argnum=argnum+2;
            case {'plot_last_batch'}
                plot_nanogrid_final=varargin{argnum+1};
                argnum=argnum+2;
        end
    end
end

B=isfield(nanogrids{numngs+1},'finDpt');
switch B
    case 1
        dpoints=nanogrids{numngs+1}.finDpt(:,1:2);
        apoints=nanogrids{numngs+1}.finApt(:,1:2);
        round2=1;
    case 0
        dpoints=nanogrids{numngs+1}.totDpt(:,1:2);
        apoints=nanogrids{numngs+1}.totApt(:,1:2);
        round2=0;
end

a=0;
for gg=1:numngs
    check=0;
    if round2==0 && ~isempty(nanogrids{gg}.dchCrd) && ~isempty(nanogrids{gg}.achCrd)
        movingtPt=nanogrids{gg}.dchCrd(:,1:2);
        refApt=nanogrids{gg}.achCrd(:,1:2);
        check=1;
    elseif round2==1 && isfield(nanogrids{gg},'redDpt') && isfield(nanogrids{gg},'redApt')...
            && ~isempty(nanogrids{gg}.redDpt) && ~isempty(nanogrids{gg}.redApt)
        movingtPt=nanogrids{gg}.redDpt(:,1:2);
        refApt=nanogrids{gg}.redApt(:,1:2);
        check=1;
    end


    if check==1
        redApt=[];
        redDpt=[];
        tmpDpt=[];
        tmpApt=[];
        mvApt=[];
        mvDpt=[];
        regDpt=[];
        regApt=[];

        %Loop through the points. For each point group surrounding points that are within 30 pixels of
        %points
        for b=1:size(movingtPt,1)
            szb=0;
            cmv=0;

            tmpD=movingtPt(b,:);
            tmpA=refApt(b,:);
            %[posind,~]= find(dpoints(:,2)>= tmpD(1,2)-rectangle & dpoints(:,2) <= tmpD(1,2)+rectangle & dpoints(:,1)>= tmpD(1,1)-rectangle & dpoints(:,1)<= tmpD(1,1)+rectangle);

            while szb<minNumPt
                [posind,~]= find(dpoints(:,2)>= tmpD(1,2)-(rectangle+cmv) & dpoints(:,2) <= tmpD(1,2)+(rectangle+cmv) & dpoints(:,1)>= tmpD(1,1)-(rectangle+cmv) & dpoints(:,1)<= tmpD(1,1)+(rectangle+cmv));
                tmpNumPt=dpoints(posind,1:2);%dchan points from reference datasets
                szb=size(tmpNumPt,1);
                cmv=cmv+1;
            end

            dpt=dpoints(posind,1:2);
            apt= apoints(posind,1:2);
            distmat=pdist2(tmpD,dpt); %calculate the pairwise distances between points in dchn and achn
            [~,ind]=sort(distmat,2);
            mappingpoints=[];
            mappingpoints(:,3:4)=dpt(ind(1,2:size(dpt,1)-1),:);
            mappingpoints(:,9:10)=apt(ind(1,2:size(dpt,1)-1),:);

            fitparmx21=polyfit(mappingpoints(:,3),mappingpoints(:,9),1);
            % Form a cell array, first member is a matrix of the
            % x1y1 pairs
            inarray{1}=[ mappingpoints(:,3) mappingpoints(:,4)];
            % second member is a vector of the output
            % x2 points
            inarray{2} = mappingpoints(:,9);
            % Input guess is [mxx21 mxy21 bx] with
            % mxy21 = 0 at first
            %****fitparmx21more=mappingfit(inarray,[fitparmx21(1) 0 fitparmx21(2) ]);
            fitparmx21more=mappingfit_fminsearch(inarray,[fitparmx21(1) 0 fitparmx21(2) ]);
            fitparmy21=polyfit(mappingpoints(:,4),mappingpoints(:,10),1);
            % Form a cell array, first member is a
            % matrix of the x1y1 pairs
            inarray{1}=[ mappingpoints(:,3) mappingpoints(:,4)];
            % second member is a vector of the output
            % y2 points
            inarray{2} = mappingpoints(:,10);
            % Input guess is [myx21 myy21 bx] with
            % myx21 = 0 at first
            %***fitparmy21more=mappingfit(inarray,[0 fitparmx21(1) fitparmx21(2) ]);
            fitparmy21more=mappingfit_fminsearch(inarray,[0 fitparmy21(1) fitparmy21(2) ]);
 
            mvdpt(:,1)=mappingfunc(fitparmx21more,tmpD);
            mvdpt(:,2)=mappingfunc(fitparmy21more,tmpD);
       
            if double(sqrt((tmpA(1,1)-mvdpt(1,1))^2 + (tmpA(1,2)-mvdpt(1,2))^2)) < minDistanceBetweenPointsA
                tmpDpt=vertcat(tmpDpt,tmpD);
                tmpApt=vertcat(tmpApt,tmpA);
            end
            if double(sqrt((tmpA(1,1)-mvdpt(1,1))^2 + (tmpA(1,2)-mvdpt(1,2))^2)) >minDistanceBetweenPointsA
                regDpt=vertcat(regDpt,tmpD);
                regApt=vertcat(regApt,tmpA);
            end
        end

        for b=1:size(tmpApt,1)
            szb=0;
            cmv=0;
            tmpA=tmpApt(b,:);
            tmpD=tmpDpt(b,:);
            % [posind,~]= find(dpoints(:,2)>= tmpA(1,2)-rectangle & dpoints(:,2) <= tmpA(1,2)+rectangle & dpoints(:,1)>= tmpA(1,1)-rectangle & dpoints(:,1)<= tmpA(1,1)+rectangle);
            while szb<minNumPt
                [posind,~]= find(dpoints(:,2)>= tmpA(1,2)-(rectangle+cmv) & dpoints(:,2) <= tmpA(1,2)+(rectangle+cmv) & dpoints(:,1)>= tmpA(1,1)-(rectangle+cmv) & dpoints(:,1)<= tmpA(1,1)+(rectangle+cmv));
                tmpNumPt=dpoints(posind,1:2);%dchan points from reference datasets
                szb=size(tmpNumPt,1);
                cmv=cmv+1;
            end
            dpt=dpoints(posind,1:2);
            apt= apoints(posind,1:2);
            distmat=pdist2(tmpA,apt); %calculate the pairwise distances between points in dchn and achn
            [~,ind]=sort(distmat,2);
            mappingpoints=[];
            mappingpoints(:,3:4)=apt(ind(1,2:size(apt,1)-1),:);
            mappingpoints(:,9:10)=dpt(ind(1,2:size(apt,1)-1),:);
            fitparmx21=polyfit(mappingpoints(:,3),mappingpoints(:,9),1);
            % Form a cell array, first member is a matrix of the
            % x1y1 pairs
            inarray{1}=[ mappingpoints(:,3) mappingpoints(:,4)];
            % second member is a vector of the output
            % x2 points
            inarray{2} = mappingpoints(:,9);
            % Input guess is [mxx21 mxy21 bx] with
            % mxy21 = 0 at first
            %****fitparmx21more=mappingfit(inarray,[fitparmx21(1) 0 fitparmx21(2) ]);
            fitparmx21more=mappingfit_fminsearch(inarray,[fitparmx21(1) 0 fitparmx21(2) ]);
            fitparmy21=polyfit(mappingpoints(:,4),mappingpoints(:,10),1);
            % Form a cell array, first member is a
            % matrix of the x1y1 pairs
            inarray{1}=[ mappingpoints(:,3) mappingpoints(:,4)];
            % second member is a vector of the output
            % y2 points
            inarray{2} = mappingpoints(:,10);
            % Input guess is [myx21 myy21 bx] with
            % myx21 = 0 at first
            %***fitparmy21more=mappingfit(inarray,[0 fitparmx21(1) fitparmx21(2) ]);
            fitparmy21more=mappingfit_fminsearch(inarray,[0 fitparmy21(1) fitparmy21(2) ]);
            mvapt(:,1)=mappingfunc(fitparmx21more,tmpA);
            mvapt(:,2)=mappingfunc(fitparmy21more,tmpA);

            if double(sqrt((tmpD(1,1)-mvapt(1,1))^2 + (tmpD(1,2)-mvapt(1,2))^2)) < minDistanceBetweenPointsA
                redDpt=vertcat(redDpt,tmpD);
                redApt=vertcat(redApt,tmpA);
                mvApt=vertcat(mvApt,mvapt);
            end
            if double(sqrt((tmpD(1,1)-mvapt(1,1))^2 + (tmpD(1,2)-mvapt(1,2))^2)) >minDistanceBetweenPointsA
                regDpt=vertcat(regDpt,tmpD);
                regApt=vertcat(regApt,tmpA);
            end


        end
        nanogrids{gg}.redDpt=redDpt;
        nanogrids{gg}.redApt=redApt;
        nanogrids{gg}.mvApt=mvApt;
        nanogrids{gg}.regApt=regApt;
        nanogrids{gg}.regDpt=regDpt;

        for b=1:size(redDpt,1)
            szb=0;
            cmv=0;
            tmpD=redDpt(b,:);
            tmpA=redApt(b,:);
            % [posind,~]= find(dpoints(:,2)>= tmpD(1,2)-rectangle & dpoints(:,2) <= tmpD(1,2)+rectangle & dpoints(:,1)>= tmpD(1,1)-rectangle & dpoints(:,1)<= tmpD(1,1)+rectangle);
            while szb<minNumPt
                [posind,~]= find(dpoints(:,2)>= tmpD(1,2)-(rectangle+cmv) & dpoints(:,2) <= tmpD(1,2)+(rectangle+cmv) & dpoints(:,1)>= tmpD(1,1)-(rectangle+cmv) & dpoints(:,1)<= tmpD(1,1)+(rectangle+cmv));
                tmpNumPt=dpoints(posind,1:2);%dchan points from reference datasets
                szb=size(tmpNumPt,1);
                cmv=cmv+1;
            end
            dpt=dpoints(posind,1:2);
            apt= apoints(posind,1:2);
            distmat=pdist2(tmpD,dpt); %calculate the pairwise distances between points in dchn and achn
            [~,ind]=sort(distmat,2);
            mappingpoints=[];
            mappingpoints(:,3:4)=dpt(ind(1,2:size(dpt,1)-1),:);
            mappingpoints(:,9:10)=apt(ind(1,2:size(dpt,1)-1),:);

            fitparmx21=polyfit(mappingpoints(:,3),mappingpoints(:,9),1);
            % Form a cell array, first member is a matrix of the
            % x1y1 pairs
            inarray{1}=[ mappingpoints(:,3) mappingpoints(:,4)];
            % second member is a vector of the output
            % x2 points
            inarray{2} = mappingpoints(:,9);
            % Input guess is [mxx21 mxy21 bx] with
            % mxy21 = 0 at first
            %****fitparmx21more=mappingfit(inarray,[fitparmx21(1) 0 fitparmx21(2) ]);
            fitparmx21more=mappingfit_fminsearch(inarray,[fitparmx21(1) 0 fitparmx21(2) ]);
            fitparmy21=polyfit(mappingpoints(:,4),mappingpoints(:,10),1);
            % Form a cell array, first member is a
            % matrix of the x1y1 pairs
            inarray{1}=[ mappingpoints(:,3) mappingpoints(:,4)];
            % second member is a vector of the output
            % y2 points
            inarray{2} = mappingpoints(:,10);
            % Input guess is [myx21 myy21 bx] with
            % myx21 = 0 at first
            %***fitparmy21more=mappingfit(inarray,[0 fitparmx21(1) fitparmx21(2) ]);
            fitparmy21more=mappingfit_fminsearch(inarray,[0 fitparmy21(1) fitparmy21(2) ]);
            mvdpt(:,1)=mappingfunc(fitparmx21more,tmpD);
            mvdpt(:,2)=mappingfunc(fitparmy21more,tmpD);
     
            if double(sqrt((tmpA(1,1)-mvdpt(1,1))^2 + (tmpA(1,2)-mvdpt(1,2))^2)) < minDistanceBetweenPointsB
                mvDpt=vertcat(mvDpt,mvdpt);
            end
        end
        nanogrids{gg}.mvDpt=mvDpt;
        if plot_final==true
            if plot_nanogrid_final==true
                if a==0

                    if ~isempty(refApt) && ~isempty(movingtPt)
                        f = figure;
                        set(f, 'Visible', 'on'); clf;
                        set(f,'Position',[624 82 896 905]);
                        ax1=subplot(1,2,1);
                        % Create scatter
                        ha1=scatter(refApt(:,1), refApt(:,2), 30, 'r'); hold on;
                        % Create scatter
                        ha2=scatter(movingtPt(:,1), movingtPt(:,2), 20, 'b','MarkerFaceColor','flat','MarkerEdgeColor','none');
                        xlim(ax1,[0 256]);
                        % Uncomment the following line to preserve the Y-limits of the axes
                        ylim(ax1,[0 512]);
                        box(ax1,'on');
                        axis(ax1,'ij');
                        % Set the remaining axes properties
                        set(ax1,'FontName','Arial Black','FontSize',16,'FontWeight','bold',...
                            'LineWidth',3,'XColor',[0 0 0],'XTick',[0 100 200],'YColor',[0 0 0],'YTick',...
                            [0 100 200 300 400 500],'ZColor',[0 0 0]);
                        hold off;
                        if ~isempty(refApt) && ~isempty(mvDpt)
                            % createfigure_overlay(refApt(:,1), refApt(:,2), 30, 'r', movingtPt(:,1), movingtPt(:,2), 20, 'b');

                            ax2=subplot(1,2,2);
                            % Create scatter
                            hb1=scatter(refApt(:,1), refApt(:,2), 30, 'r'); hold on;
                            % Create scatter
                            hb2=scatter(mvDpt(:,1), mvDpt(:,2), 20, 'g','MarkerFaceColor','flat','MarkerEdgeColor','none');
                            xlim(ax2,[0 256]);
                            % Uncomment the following line to preserve the Y-limits of the axes
                            ylim(ax2,[0 512]);
                            box(ax2,'on');
                            axis(ax2,'ij');
                            % Set the remaining axes properties
                            set(ax2,'FontName','Arial Black','FontSize',16,'FontWeight','bold',...
                                'LineWidth',3,'XColor',[0 0 0],'XTick',[0 100 200],'YColor',[0 0 0],'YTick',...
                                [0 100 200 300 400 500],'ZColor',[0 0 0]);
                            movegui('center')
                            saveas(gcf, ['Points_overLap_' num2str(gg) '.png']);
                            saveas(gcf, ['Points_overLap_' num2str(gg) '.svg']);
                            hold off;
                            a=1;

                        else
                            ax2=subplot(1,2,2);
                            % Create scatter
                            hb1=scatter(1, 1, 30, 'r'); hold on;
                            % Create scatter
                            hb2=scatter(1, 1, 20, 'g','MarkerFaceColor','flat','MarkerEdgeColor','none');
                            xlim(ax2,[0 256]);
                            % Uncomment the following line to preserve the Y-limits of the axes
                            ylim(ax2,[0 512]);
                            box(ax2,'on');
                            axis(ax2,'ij');
                            % Set the remaining axes properties
                            set(ax2,'FontName','Arial Black','FontSize',16,'FontWeight','bold',...
                                'LineWidth',3,'XColor',[0 0 0],'XTick',[0 100 200],'YColor',[0 0 0],'YTick',...
                                [0 100 200 300 400 500],'ZColor',[0 0 0]);
                            movegui('center')
                            saveas(gcf, ['Points_overLap_' num2str(gg) '.png']);
                            saveas(gcf, ['Points_overLap_' num2str(gg) '.svg']);
                            hold off;
                            a=1;



                        end
                    end
                end
                if a==1
                    if ~isempty(refApt) && ~isempty(movingtPt)
                        set(ha1,'XData',refApt(:,1),'YData',refApt(:,2));
                        set(ha2,'XData',movingtPt(:,1),'YData',movingtPt(:,2));
                        if ~isempty(refApt) && ~isempty(mvDpt)
                            set(hb1,'XData',refApt(:,1),'YData',refApt(:,2));
                            set(hb2,'XData',mvDpt(:,1),'YData',mvDpt(:,2));
                            saveas(gcf, ['Points_overLap_' num2str(gg) '.png']);
                            saveas(gcf, ['Points_overLap_' num2str(gg) '.svg']);

                        else
                            set(hb1,'XData',1,'YData',1);
                            set(hb2,'XData',1,'YData',1);
                            saveas(gcf, ['Points_overLap_' num2str(gg) '.png']);
                            saveas(gcf, ['Points_overLap_' num2str(gg) '.svg']);
                        end

                    end


                end
            end
        end
        if plot_each_image==true
            if a==0

                if ~isempty(refApt) && ~isempty(movingtPt)
                    % createfigure_overlay(refApt(:,1), refApt(:,2), 30, 'r', mvDpt(:,1), mvDpt(:,2), 20, 'g');
                    %createfigure_overlay(X1, Y1, S1, C1, X2, Y2, S2, C2)
                    % if plot_each_image==true
                    f = figure;
                    set(f, 'Visible', 'on'); clf;
                    set(f,'Position',[624 82 896 905]);
                    ax1=subplot(1,2,1);
                    % Create scatter
                    ha1=scatter(refApt(:,1), refApt(:,2), 30, 'r'); hold on;
                    % Create scatter
                    ha2=scatter(movingtPt(:,1), movingtPt(:,2), 20, 'b','MarkerFaceColor','flat','MarkerEdgeColor','none');
                    xlim(ax1,[0 256]);
                    % Uncomment the following line to preserve the Y-limits of the axes
                    ylim(ax1,[0 512]);
                    box(ax1,'on');
                    axis(ax1,'ij');
                    % Set the remaining axes properties
                    set(ax1,'FontName','Arial Black','FontSize',16,'FontWeight','bold',...
                        'LineWidth',3,'XColor',[0 0 0],'XTick',[0 100 200],'YColor',[0 0 0],'YTick',...
                        [0 100 200 300 400 500],'ZColor',[0 0 0]);
                    hold off;
                    if ~isempty(refApt) && ~isempty(mvDpt)
                        % createfigure_overlay(refApt(:,1), refApt(:,2), 30, 'r', movingtPt(:,1), movingtPt(:,2), 20, 'b');

                        ax2=subplot(1,2,2);
                        % Create scatter
                        hb1=scatter(refApt(:,1), refApt(:,2), 30, 'r'); hold on;
                        % Create scatter
                        hb2=scatter(mvDpt(:,1), mvDpt(:,2), 20, 'g','MarkerFaceColor','flat','MarkerEdgeColor','none');
                        xlim(ax2,[0 256]);
                        % Uncomment the following line to preserve the Y-limits of the axes
                        ylim(ax2,[0 512]);
                        box(ax2,'on');
                        axis(ax2,'ij');
                        % Set the remaining axes properties
                        set(ax2,'FontName','Arial Black','FontSize',16,'FontWeight','bold',...
                            'LineWidth',3,'XColor',[0 0 0],'XTick',[0 100 200],'YColor',[0 0 0],'YTick',...
                            [0 100 200 300 400 500],'ZColor',[0 0 0]);
                        movegui('center')
                        saveas(gcf, ['Points_overLap_' num2str(gg) '.png']);
                        hold off;
                        a=1;

                    else
                        ax2=subplot(1,2,2);
                        % Create scatter
                        hb1=scatter(1, 1, 30, 'r'); hold on;
                        % Create scatter
                        hb2=scatter(1, 1, 20, 'g','MarkerFaceColor','flat','MarkerEdgeColor','none');
                        xlim(ax2,[0 256]);
                        % Uncomment the following line to preserve the Y-limits of the axes
                        ylim(ax2,[0 512]);
                        box(ax2,'on');
                        axis(ax2,'ij');
                        % Set the remaining axes properties
                        set(ax2,'FontName','Arial Black','FontSize',16,'FontWeight','bold',...
                            'LineWidth',3,'XColor',[0 0 0],'XTick',[0 100 200],'YColor',[0 0 0],'YTick',...
                            [0 100 200 300 400 500],'ZColor',[0 0 0]);
                        movegui('center')
                        saveas(gcf, ['Points_overLap_' num2str(gg) '.png']);
                        hold off;
                        a=1;



                    end
                end
            end
            if a==1
                if ~isempty(refApt) && ~isempty(movingtPt)
                    set(ha1,'XData',refApt(:,1),'YData',refApt(:,2));
                    set(ha2,'XData',movingtPt(:,1),'YData',movingtPt(:,2));
                    if ~isempty(refApt) && ~isempty(mvDpt)
                        set(hb1,'XData',refApt(:,1),'YData',refApt(:,2));
                        set(hb2,'XData',mvDpt(:,1),'YData',mvDpt(:,2));
                        saveas(gcf, ['Points_overLap_' num2str(gg) '.png']);

                    else
                        set(hb1,'XData',1,'YData',1);
                        set(hb2,'XData',1,'YData',1);
                        saveas(gcf, ['Points_overLap_' num2str(gg) '.png']);

                    end

                end


            end
        end



    end

end


finDpt=[];
finApt=[];
mvallApt=[];
mvallDpt=[];
regallApt=[];
regallDpt=[];
for x=1:numngs
    if ~isempty(nanogrids{x}.dchCrd) && ~isempty(nanogrids{x}.achCrd)

        finDpt=vertcat(finDpt,nanogrids{x}.redDpt);
        finApt=vertcat(finApt,nanogrids{x}.redApt);
        mvallApt=vertcat(mvallApt,nanogrids{x}.mvApt);
        mvallDpt=vertcat(mvallDpt,nanogrids{x}.mvDpt);
        regallDpt=vertcat(regallDpt,nanogrids{x}.regDpt);
        regallApt=vertcat(regallApt,nanogrids{x}.regApt);
    end
end
if plot_final==true
    createfigure_overlay(finApt(:,1), finApt(:,2), 10, 'r', mvallDpt(:,1), mvallDpt(:,2), 5, 'g');
    saveas(gcf, ['Points_overLap_All_points_cutoff_' num2str(CutoffDistance) '_pixels.png']);
     saveas(gcf, ['Points_overLap_All_points_cutoff_' num2str(CutoffDistance) '_pixels.svg']);
    createfigure_overlay(apoints(:,1), apoints(:,2), 10, 'r', dpoints(:,1), dpoints(:,2), 5, 'b');
    saveas(gcf, ['Points_Initial_All_points_cutoff_' num2str(CutoffDistance) '_pixels.png']);
     saveas(gcf, ['Points_Initial_All_points_cutoff_' num2str(CutoffDistance) '_pixels.svg']);
    if ~isempty(regallApt) & ~isempty(regallDpt)
        createfigure_overlay(regallApt(:,1), regallApt(:,2), 30, 'r', regallDpt(:,1), regallDpt(:,2), 20, 'b');
        saveas(gcf, ['Regect_Initial_All_points_cutoff_' num2str(CutoffDistance) '_pixels.png']);
        saveas(gcf, ['Regect_Initial_All_points_cutoff_' num2str(CutoffDistance) '_pixels.svg']);
    end
end
nanogrids{numngs+1}.finDpt=finDpt;
nanogrids{numngs+1}.finApt=finApt;
nanogrids{numngs+1}.mvallApt=mvallApt;
nanogrids{numngs+1}.mvallDpt=mvallDpt;
nanogrids{numngs+1}.regallApt=regallApt;
nanogrids{numngs+1}.regallDpt=regallDpt;
finalRef.finDpt=finDpt;
finalRef.finApt=finApt;
close all;
filepath=fullfile(path_to_data,[outf]);
save([outf '_finalRef'],'finalRef')
save(filepath,'finalRef')
save(outf,'nanogrids')
end