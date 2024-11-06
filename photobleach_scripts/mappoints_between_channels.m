function [transpoints]=mappoints_between_channels(pointsToMap,refPointsA,refPointsB,rectangle,number_Points_to_Fit) 
%[transpoints]=mappoints(pointsToMap,refPointsA,refPointsB) 
%This funtion will map points (pointsToMap) from channelA to channelB based
%on Transformation between reference points in channelA (refPointsA) and
%reference points in channelB (refPointsB). To improve the quality of the
%fit, this function uses reference points that are within a window (rectangle) of movingpoints   
        ArefPt=refPointsA; 
        BrefPts=refPointsB;
        transpt=[];
        movingtPt=pointsToMap;
        minNumPt=number_Points_to_Fit;
        %Loop through the points. For each point group surrounding points that are within 30 pixels of
        %points
                    for b=1:size(movingtPt,1)
                        cmv=0;
                        szb=0;
                        tmp=movingtPt(b,:);
                        %tmpb=refApt(b,:);
                        while szb<minNumPt
                            [posind,~]= find(ArefPt(:,2)>= tmp(1,2)-(rectangle+cmv) & ArefPt(:,2) <= tmp(1,2)+(rectangle+cmv) & ArefPt(:,1)>= tmp(1,1)-(rectangle+cmv) & ArefPt(:,1)<= tmp(1,1)+(rectangle+cmv));
                            tmpNumPt=ArefPt(posind,1:2);%dchan points from reference datasets
                            szb=size(tmpNumPt,1);
                            cmv=cmv+1;
                        end
                        AwndwPt=ArefPt(posind,1:2);%Cy3 points from reference datasets
                        BwndwPt= BrefPts(posind,1:2);%Cy5 points from reference datasets
                        distmat=pdist2(tmp,AwndwPt); %calculate the pairwise distances between points in dchn and achn
                        [~,ind]=sort(distmat,2);
                        mappingpoints=[];
                        mappingpoints(:,3:4)=AwndwPt(ind(1,2:size(AwndwPt,1)-1),:);
                        mappingpoints(:,9:10)=BwndwPt(ind(1,2:size(AwndwPt,1)-1),:);

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

                        fitparmy21more=mappingfit_fminsearch(inarray,[0 fitparmy21(1) fitparmy21(2) ]);


                        mvdpt(:,1)=mappingfunc(fitparmx21more,tmp);
                        mvdpt(:,2)=mappingfunc(fitparmy21more,tmp);
                        transpt=vertcat(transpt,mvdpt);
                    end
                    transpoints=transpt;
end