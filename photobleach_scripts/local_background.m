function [bkgImg, stdImg] = local_background(Image, HalfBox, Blocksearch, Border)
    % local_background computes the local background and standard deviation images.
    %
    % [bkgImg, stdImg] = local_background(Image, HalfBox, Blocksearch, Border)
    %
    % This function generates two images:
    % - bkgImg: local background image (median of pixel intensities in blocks)
    % - stdImg: local standard deviation image (standard deviation of pixel intensities in blocks)
    %
    % Parameters:
    % Image     : The input image for local background computation.
    % HalfBox   : The step size to move the window across the image.
    % Blocksearch : The size of the block to compute statistics over.
    % Border    : The width of the image border to exclude from processing.
    %
    % Example:
    % [bkgImg, stdImg] = local_background(myImage, 10, 20, 5)
    %
    % Copyright 2024 Andrew Moreno, Harvard Medical School.
    % Licensed under the GNU General Public License, version 3.

    % Assign input parameters to local variables for clarity
    img = Image;
    Stepper = HalfBox;
    Blocksize = Blocksearch;

    % Define the range of rows and columns to process, excluding the border
    Rowcent = Border:Stepper:(512 - Border);
    Colcent = Border:Stepper:(256 - Border);

    % Initialize the output images
    bkgimg = zeros(512, 256);
    bkgimgSTD = zeros(512, 256);

    % Iterate over the defined rows and columns in steps of 2
    for jj = 1:2:size(Rowcent, 2)
        for ff = 1:2:size(Colcent, 2)
            % Define the limits of the block to be processed
            [BXxlow, BXxhi, BXylow, BXyhi] = ROI_boundry([Rowcent(jj), Colcent(ff)], Blocksize / 2);
            BXwndw = double(img(BXxlow:BXxhi, BXylow:BXyhi));
            block_intensity = BXwndw(:); % Create a vector of the block's intensity values

            % Outlier removal process
            zFactor = 1; % Threshold for outlier detection
            outliers = true; % Initialize outliers flag

            while any(outliers)
                % Compute standard deviation and median of the block
                stdDev = std(block_intensity);
                medianValue = median(block_intensity);

                % Identify outliers based on the zFactor threshold
                outliers = abs(block_intensity - medianValue) > (zFactor * stdDev);

                % Remove outliers from the block intensity vector
                block_intensity = block_intensity(~outliers);
            end

            % Define the limits of the smaller block for assigning the median and std deviation
            [xlow, xhi, ylow, yhi] = ROI_boundry([Rowcent(jj), Colcent(ff)], Stepper);

            % Assign the computed median and standard deviation to the output images
            bkgimg(xlow:xhi, ylow:yhi) = medianValue;
            bkgimgSTD(xlow:xhi, ylow:yhi) = stdDev;
        end
    end

    % Return the computed background and standard deviation images
    bkgImg = bkgimg;
    stdImg = bkgimgSTD;

end

function [xlow, xhi, ylow, yhi] = ROI_boundry(center, halfSize)
% ROI_boundry computes the limits of a block centered at 'center' with a size of 'halfSize'.
%
% [xlow, xhi, ylow, yhi] = ROI_boundry(center, halfSize)
%
% center  : [x, y] coordinates of the center of the block.
% halfSize: Half the size of the block.
%
% Returns:
% xlow, xhi: The lower and upper limits in the x-direction.
% ylow, yhi: The lower and upper limits in the y-direction.

xlow=round(center(1)-halfSize+0.5);
xhi=round(xlow+2*halfSize-1);
ylow=round(center(2)-halfSize+0.5);
yhi=round(ylow+2*halfSize-1);
end
