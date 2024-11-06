function make_boxplot(Goodpoints,Badpoints,boxsize,text)

rjpts=Badpoints;
kppts=Goodpoints;
spotbx=boxsize;
yestxt=text;
        
for gg=1:size(rjpts,1)
    D=@(i) box_point(rjpts(i,1:2),(spotbx-1)/2,(spotbx-1)/2,'r');
    D(gg);
    if yestxt==1
        F=@(h) text(rjpts(h,1),rjpts(h,2),num2str(h),'Color',[0.7 0.2 0.7]);
        F(gg);
    end
    hold on;
end
for gg=1:size(kppts,1)
    JJ=@(i) box_point(kppts(i,1:2),(spotbx-1)/2,(spotbx-1)/2,'g');
    JJ(gg);
    if yestxt==1
        F=@(h) text(kppts(h,1),kppts(h,2),num2str(h),'Color',[0.7 0.2 .7]);
        F(gg);
    end
    hold on;
end
set(gca,'FontName','Arial Black','FontSize',16,'FontWeight','bold','LineWidth',3)%,'XColor',[0 0 0],'YColor',[0 0 0],'YTick','ZColor',[0 0 0]);
end


function box_point(center,xdis,ydis,linecolor)
xlin = [center(1)-xdis  center(1)+xdis  center(1)+xdis  center(1)-xdis center(1)-xdis];
ylin =[center(2)-ydis  center(2)-ydis  center(2)+ydis  center(2)+ydis center(2)-ydis] ;
line(xlin,ylin,'color',linecolor,'LineWidth',1)
end