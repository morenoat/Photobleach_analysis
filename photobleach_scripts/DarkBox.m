function [fnpt]=DarkBox(image,border,minDist,brightspots,boxwidth,FinalPoints,segmentnumber,outname)
%%Function negativespace(image,border,window)
% This function will identify areas (window) of an image (image) that lack a
% fluorescent molecule. The code will then compile a list of
% corridinates (xy) corresponding to the center of the windows that do
% not "touch" each other
    dataout=[pwd '/' outname];
    aIout=image;
    alldpts=brightspots(brightspots(:,8)>0,:);
%Creating an array that will specify the start and stop of each grid 
se=strel('disk',boxwidth,0);
sen=se.getnhood();
senReal=sen.*1;
senReal(boxwidth+1,boxwidth+1)=9;
BGbx=sen-1;
BGbx(BGbx==-1)=1;
 IntNum=sum(sen(:));

window=size(senReal,1);
M=ceil(512/window);
N=ceil(256/window);
newMap=repmat(senReal,M,N);
mrklocalA=[];
for jk=1:size(newMap,1)
    spt=[];
    spt=(find(newMap(jk,:)==9))';
    if~isempty(spt)
        spt(:,2)=jk;
        mrklocalA=vertcat(mrklocalA,spt);
    end
end
posind=mrklocalA;
posind=posind(posind(:,1)>border & posind(:,2)>border & posind(:,2)<512-border & posind(:,1)<256-border,:);%keep boxes that are within border
    gdbx(:,1)=posind(:,1)-boxwidth;
    gdbx(:,2)=posind(:,2)-boxwidth;
    gdbx(:,3)=posind(:,1)+boxwidth;
    gdbx(:,4)=posind(:,2)+boxwidth;
    cnpt=posind;

    dDistmat=pdist2(cnpt(:,1:2),alldpts(:,1:2)); %calculate the pairwise distances between points in darkspots and dchn brightspots 
    [blkval,~]=sort(dDistmat,2); %sort the distances down columns 
    gdbxB=gdbx(blkval(:,1)>minDist,:); %remove spots that are inside distance threshold leading to 'spurious detection'
    cnptB=cnpt(blkval(:,1)>minDist,:); %remove spots that are inside distance threshold leading to 'spurious detection'

%% calculate the integrated intensity for background aareas
    pts=cnptB;
    tmpBG=[];
   for ff=1:size(pts,1)
        pc2=zeros(1,6);
        pc2(1,1)=pts(ff,1);
        pc2(1,2)=pts(ff,2);
        aoix=pts(ff,2);
        aoiy=pts(ff,1);
        [xlow,xhi,ylow,yhi]=AOI_Limits([aoix aoiy],window/2);
        wndw=aIout(xlow:xhi,ylow:yhi);
        box = wndw.*sen;
        boxBg=wndw.*BGbx;
        ptI = sum(box(:));
        bgI = median((boxBg(boxBg(:)~=0)));%median background intensity
        pc2(1,3)=ptI;
        pc2(1,4)=bgI;
        pc2(1,5)=bgI*(IntNum);
        pc2(1,6)=ptI-(bgI*(IntNum)); 
        tmpBG(ff,:)=pc2'; 
   end
    [fxIpts(:,2),fxIpts(:,1)]=ecdf(alldpts(:,8));
    [~ ,idx]=min(abs(fxIpts(:,2)-.0075));
    cutI=fxIpts(idx,1);
    bdpt=[];
    bdpt=gdbxB(tmpBG(:,6)>cutI,:);
    cnptB(tmpBG(:,6)>cutI,:)=[];
    gdbxB(tmpBG(:,6)>cutI,:)=[];
    



    
%keep boxes such that none are touching 
    tmp=[0 0];
    fnbx=[];
    
    fnpt=[];
    idGood=[];
    idExld=[];
    cBD=1;
    cGD=1;
    tmpbdpt=[];

    if size(cnptB,1)<size(FinalPoints,1)
        fnbx=gdbxB;
        fnpt=cnptB;

    else
        for dd=1:size(cnptB,1)
            if ~any(ismember(tmp,cnptB(dd,:),'rows'))
                fnbx=vertcat(fnbx,gdbxB(dd,:));
                fnpt=vertcat(fnpt,cnptB(dd,:));
                distmat=pdist2(cnptB(dd,:),cnptB);
                [val,ind]=sort(distmat,2);
                tmp=vertcat(tmp,cnptB(ind(:,val(1,:)==window),:));
                idGood(cGD,1)=dd;
                cGD=cGD+1;
            else
                idExld(cBD,1)=dd;
                cBD=cBD+1;
            end
        end
        if size(idGood,1)<size(FinalPoints,1)
            numdiff=size(FinalPoints,1)-size(idGood,1);
            addtoLst=(randperm(size(idExld,1),numdiff))';
            idAdd=idExld(addtoLst,1);
            idExld(addtoLst,:)=[];
            AddCrd=gdbxB(idAdd,:);
            AddPTS=cnpt(idAdd,:);
            fnbx=vertcat(fnbx,AddCrd);
            fnpt=vertcat(fnpt,AddPTS);
            tmpbdpt=gdbxB(idExld,:);

        else
             tmpbdpt=gdbxB(idExld,:);

        end
    end
    bdpt=vertcat(bdpt,tmpbdpt);

    figure;
    set(gcf,'position',[696 29 637 952]);
    imagesc(image);
    hold on;
%     for gg=1:size(gdbx,1)
%         line([gdbx(gg,2) gdbx(gg,4)],[gdbx(gg,1) gdbx(gg,1)],'LineWidth',1,'Color','red');
%         line([gdbx(gg,2) gdbx(gg,4)],[gdbx(gg,3) gdbx(gg,3)],'LineWidth',1,'Color','red');
% 
%         line([gdbx(gg,2) gdbx(gg,2)],[gdbx(gg,1) gdbx(gg,3)],'LineWidth',1,'Color','red');
%         line([gdbx(gg,4) gdbx(gg,4)],[gdbx(gg,1) gdbx(gg,3)],'LineWidth',1,'Color','red');    
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Make boxes of windows green=final negative boxes red=rejected boxes%%%%%%%%%%%%%
    for gg=1:size(bdpt,1)
        line([bdpt(gg,1) bdpt(gg,1)],[bdpt(gg,2) bdpt(gg,4)],'LineWidth',2,'Color','r');
        line([bdpt(gg,3) bdpt(gg,3)],[bdpt(gg,2) bdpt(gg,4)],'LineWidth',2,'Color','r');

        line([bdpt(gg,1) bdpt(gg,3)],[bdpt(gg,2) bdpt(gg,2)],'LineWidth',2,'Color','r');
        line([bdpt(gg,1) bdpt(gg,3)],[bdpt(gg,4) bdpt(gg,4)],'LineWidth',2,'Color','r');    
    end
    for gg=1:size(fnbx,1)
        line([fnbx(gg,1) fnbx(gg,1)],[fnbx(gg,2) fnbx(gg,4)],'LineWidth',2,'Color','g');
        line([fnbx(gg,3) fnbx(gg,3)],[fnbx(gg,2) fnbx(gg,4)],'LineWidth',2,'Color','g');

        line([fnbx(gg,1) fnbx(gg,3)],[fnbx(gg,2) fnbx(gg,2)],'LineWidth',2,'Color','g');
        line([fnbx(gg,1) fnbx(gg,3)],[fnbx(gg,4) fnbx(gg,4)],'LineWidth',2,'Color','g');    
    end
    movegui center

     set(gca,'FontName','Arial Black','FontSize',16,'FontWeight','bold','LineWidth',3)%,'XColor',[0 0 0],'YColor',[0 0 0],'YTick','ZColor',[0 0 0]);
     saveas(gcf, fullfile(dataout,[outname num2str(segmentnumber) '_Negativespace']), 'png');hold off;
     close all;
end
  