function [bestMol,removeMol,MedVal,StdVal]=curateList_simple(molecules,ColToCur,Deviate)


       
zFactor = Deviate; % or whatever you want
outliers=1;
A=molecules(~isnan(molecules(:,ColToCur)),:);
Outlist=molecules(isnan(molecules(:,ColToCur)),:);%molecules that will be removed based on deviating from median
while any(outliers)==1
    tmp=[];
    outliers=[];
    stdDev = std(A(:,ColToCur)); % Compute standard deviation
    medianValue = median(A(:,ColToCur)); % Compute median
    %Create a binary map of where outliers live
    outliers = abs(A(:,ColToCur)-medianValue) > (zFactor * stdDev);
    tmp=A(outliers==1,:);
    Outlist=vertcat(Outlist,tmp);
    A=A(outliers==0,:);
    
    stdDev = std(A(:,ColToCur)); % Compute standard deviation
    medianValue = median(A(:,ColToCur)); % Compute median
    %Create a binary map of where outliers live
    outliers = abs(A(:,ColToCur)-medianValue) > (zFactor * stdDev);
end
bestMol=A;
removeMol=Outlist;
MedVal=medianValue;
StdVal=stdDev;
end