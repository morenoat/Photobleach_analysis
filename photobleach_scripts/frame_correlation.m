function [fcorrmatLeftSide,fcorrmatRightSide] = frame_correlation(reader,period,firstCy3,border)


% Uses only the left half of the image (Cy3 channel) to calculate image 
% correlation
%Also provides right side image correlation this is output as correlation
%over periods starting at firstCy3 image
%

% nframes = reader.getSizeT;
nframes = reader.getImageCount;

fcorrmatLeftSide = [];
fcorrmatRightSide = [];
imcurr = zeros(512);

for n=firstCy3:1:(firstCy3+period-1)
    imcurr = imcurr + double(bfGetPlane(reader,n));
end

imcurr1 = imcurr(:,1+border:256-border);
 imcurr2 = imcurr(:,257+border:end-border);

imcurr_minus_mean1 = imcurr1(:) - mean(imcurr1(:));
imcurr_minus_mean2 = imcurr2(:) - mean(imcurr2(:));
imcurr_std1 = std(imcurr1(:));
imcurr_std2 = std(imcurr2(:));

w = waitbar(0,'Identifying distinct fields of view.');

for j =firstCy3+ period:period:nframes-(firstCy3+period)+1  
    
    waitbar(j/nframes,w,'Identifying distinct fields of view.');
    
    % Add up all frames in one period of the 
    imnext = double(bfGetPlane(reader,j));
    for n=1:period-1
        imnext = imnext + double(bfGetPlane(reader,j+n));
    end

    imnext1 = imnext(:,1+border:256-border);
   imnext2 = imnext(:,257+border:end-border);
    
    imnext_minus_mean1 = imnext1(:)-mean(imnext1(:));
    imnext_minus_mean2 = imnext2(:)-mean(imnext2(:));
    imnext_std1 = std(imnext1(:));
    imnext_std2 = std(imnext2(:));
    
    pxcorr1 = imnext_minus_mean1.*imcurr_minus_mean1/(imnext_std1*imcurr_std1);
    pxcorr2 = imnext_minus_mean2.*imcurr_minus_mean2/(imnext_std2*imcurr_std2);
    fcorrmatLeftSide(end+1) = mean(pxcorr1);
    fcorrmatRightSide(end+1) = mean(pxcorr2);
    
    imcurr_minus_mean1 = imnext_minus_mean1;
    imcurr_minus_mean2 = imnext_minus_mean2;
    imcurr_std1 = imnext_std1;
    imcurr_std2 = imnext_std2;
    

end

close(w);

end