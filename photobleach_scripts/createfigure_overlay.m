function createfigure_overlay(X1, Y1, S1, C1, X2, Y2, S2, C2)
%CREATEFIGURE(X1, Y1, S1, C1, X2, Y2, S2, C2)
%  X1:  scatter x
%  Y1:  scatter y
%  S1:  scatter s
%  C1:  scatter c
%  X2:  scatter x
%  Y2:  scatter y
%  S2:  scatter s
%  C2:  scatter c

%  Auto-generated by MATLAB on 26-Apr-2018 15:16:36

% Create figure
f = figure;
set(f, 'Visible', 'on'); clf;
set(f,'Position',[624 81 650 906]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create scatter
scatter(X1,Y1,S1,C1,'MarkerFaceColor','flat','MarkerEdgeColor','k');

% Create scatter
scatter(X2,Y2,S2,C2,'MarkerFaceColor','flat','MarkerEdgeColor','none');

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 256]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0 512]);
box(axes1,'on');
axis(axes1,'ij');
% Set the remaining axes properties
set(axes1,'FontName','Arial Black','FontSize',16,'FontWeight','bold',...
    'LineWidth',3,'XColor',[0 0 0],'XTick',[0 100 200],'YColor',[0 0 0],'YTick',...
    [0 100 200 300 400 500],'ZColor',[0 0 0]);
end