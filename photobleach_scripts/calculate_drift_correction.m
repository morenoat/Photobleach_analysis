function [xval,yval,fit_output]=calculate_drift_correction(Ptlist,reader,framestart,frameend,timeMovie,thres,border,LargeBox,area,SpotCut,SigRatio,Red_Green,segment,outfname,varargin)

    dataout=[pwd '\' outfname];
    limdrft=10;
    threshold=thres;
    trajlength=frameend - framestart;
    fitstart=1;
    numCy3img=[];
    c=1;
    for jj=framestart:frameend
        if timeMovie(jj,2)==1
    numCy3img(c,:)=timeMovie(jj,:);
    c=c+1;
        end
    end
    datrange=size(numCy3img,1);
    masterPt=Ptlist(:,1:2);
    inlength=0;
    red=Red_Green;
    % parameters, if present
    if inlength>0
        SG_Smooth=varargin{1}(:);                  
        SG_PolyOrderX=SG_Smooth(1);                  %
        SG_FrameX=SG_Smooth(2);                      %
        SG_PolyOrderY=SG_Smooth(3);                  %
        SG_FrameY=SG_Smooth(4);                      %                                              

    else
        SG_PolyOrderX=5;                  %
        SG_FrameX=17;                      %
        SG_PolyOrderY=5;                  %
        SG_FrameY=17;   
    end

    n=1;
    c=1;
    newmrk=[];
    xChange=[];
    yChange=[];
    for j =framestart:frameend
        if timeMovie(j,2)==1 % if this is an Cy3 frame
            newpos=[];
            oldpos=[];
            oldposB=[];
            xdiff=0;
            ydiff=0;
            doutcent=[0,0];
            currframe = double(bfGetPlane(reader,j));
            if red==1
            Iraw = currframe(:,257:end); % 641nm half of the image
            else
                Iraw = currframe(:,1:256); % 641nm half of the image
            end
            halfBox=area/2;
            MaskSpot=SpotCut;
            [bkgImg,~]=Local_Bkgrd(Iraw,halfBox,LargeBox,border);
            bkgImg(bkgImg==0)=median(bkgImg(:));
            Isub=Iraw-bkgImg;%subtract local background
            Ifit=relnoise(Isub,7,3,'disk');%local background subtraction
            Allval=Ifit(:);
            [bestbind,~,~,~]=curateList_simple(Allval,1,4);
            bgThresh=prctile(bestbind,threshold);
            %Using the FastPeakFind to determine xy coordinates of peaks in dchan and achan
            dpkRaw=double(peak_find_fast(Ifit,bgThresh));%FastPeakFind(d, thres, filt ,edg, res, fid)
            dpkRaw=[dpkRaw(1:2:end),dpkRaw(2:2:end)]; %xy coordniate of donor channel
            doutcentB=cntrd(Ifit,round(dpkRaw),5);%calculates the centroid of bright spots to sub-pixel accuracy.
            if ~isempty(doutcentB)
                Psfpt=[];
                doutcent=vertcat(doutcent,doutcentB(:,1:2));
                [mpidx,~]=find(doutcent(:,2)+area<512 & doutcent(:,2)-area>0 & doutcent(:,1)+area<256 & doutcent(:,1)-area>0);
                if ~isempty(mpidx)
                    kpInpt=doutcent(mpidx,1:2);
                    [~,Psfpt]=fitSpot_drift(Iraw,Ifit,bkgImg,kpInpt(:,1:2),area,MaskSpot,SigRatio);%Model spots as 2D gaussian for the shape and intensity of the spots
                    distmat=pdist2(masterPt,Psfpt(:,1:2));
                    [val,ind]=sort(distmat,2);
                    newpos=Psfpt(ind(val(:,1)<limdrft,1),1:2);
                    distmat=pdist2(newpos,masterPt(:,1:2));
                    [valb,indb]=sort(distmat,2);
                    oldpos=masterPt(indb(valb(:,1)<limdrft,1),1:2);
                    if ~isempty(newmrk)

                        distmat=pdist2(newmrk,Psfpt(:,1:2));
                        [val,ind]=sort(distmat,2);
                        newposB=Psfpt(ind(val(:,1)<limdrft,1),1:2);
                        distmat=pdist2(newposB,newmrk(:,1:2));
                        [valb,indb]=sort(distmat,2);
                        oldposB=newmrk(indb(valb(:,1)<limdrft,1),1:2);


                    end
                end
            end
            
            if ~isempty(oldpos) && ~isempty(newpos)
                xdiff=(oldpos(:,1)-newpos(:,1));
                ydiff=(oldpos(:,2)-newpos(:,2));
            end
            
            if ~isempty(oldposB) && ~isempty(newposB)
                xdiffB=(oldposB(:,1)-newposB(:,1));
                ydiffB=(oldposB(:,2)-newposB(:,2));
                xdiff=vertcat(xdiff,xdiffB);
                ydiff=vertcat(ydiff,ydiffB);
            end
            
            
            if ~isempty(Psfpt)
                newmrk=Psfpt;
            end
            xChange(c,1)=median(xdiff);
            yChange(c,1)=median(ydiff);

            n=n+1;
            c=c+1;
        end
    end

    valx=sgolayfilt(xChange,SG_PolyOrderX,SG_FrameX);
    valy=sgolayfilt(yChange,SG_PolyOrderY,SG_FrameY);

    f = figure;
    set(f, 'Visible', 'on'); clf;
    plot(numCy3img(:,1),xChange,'b',numCy3img(:,1),valx,'r','Linewidth',3);
    title(['xdrift, polyfit order:' num2str(SG_PolyOrderX)],'FontName','Arial Black','FontSize',16,'FontWeight','bold');
    xlabel('frame number');ylabel('x pixel')
    set(gca,'FontName','Arial Black','FontSize',16,'FontWeight','bold','LineWidth',3)%,'XColor',[0 0 0],'YColor',[0 0 0],'YTick','ZColor',[0 0 0]);
    saveas(gcf, fullfile(dataout, [outfname '_' num2str(segment) 'Cy5Xdrift_grid']), 'png');hold off; 
    saveas(gcf, fullfile(dataout, [outfname '_' num2str(segment) 'Cy5Xdrift_grid']), 'svg');hold off; 
    f = figure;
    set(f, 'Visible', 'on'); clf;
    plot(numCy3img(:,1),yChange,'b',numCy3img(:,1),valy,'r','Linewidth',3);
    title(['ydrift, polyfit order:' num2str(SG_PolyOrderY)],'FontName','Arial Black','FontSize',16,'FontWeight','bold');
    xlabel('frame number');ylabel('y pixel')
    set(gca,'FontName','Arial Black','FontSize',16,'FontWeight','bold','LineWidth',3)%,'XColor',[0 0 0],'YColor',[0 0 0],'YTick','ZColor',[0 0 0]);
    saveas(gcf, fullfile(dataout, [outfname '_' num2str(segment) 'Cy5Ydrift_grid']), 'png');hold off; 
    saveas(gcf, fullfile(dataout, [outfname '_' num2str(segment) 'Cy5Ydrift_grid']), 'svg');hold off; 
    close all;
    tmpxval=zeros(frameend,1);
    tmpyval=zeros(frameend,1);
    ggg=1;
    fit_data=[];
    fit_data=numCy3img(:,1);
    fit_data(:,2)=xChange;
    fit_data(:,3)=valx;
    fit_data(:,4)=yChange;
    fit_data(:,5)=valy;
    for j =framestart:frameend
        if timeMovie(j,2)==1
            tmpxval(j,1)=valx(ggg,1);
            tmpyval(j,1)=valy(ggg,1);
            ggg=ggg+1;
        end
        if timeMovie(j,2)~=1 && ggg>=2
            tmpxval(j,1)=valx(ggg-1,1);
            tmpyval(j,1)=valy(ggg-1,1);
        elseif timeMovie(j,2)~=1 && ggg<2
            tmpxval(j,1)=valx(ggg,1);
            tmpyval(j,1)=valy(ggg,1);
        end        
        
    end
            
    xval=tmpxval;
    yval=tmpyval;
     fit_output=fit_data;  
        
    end
