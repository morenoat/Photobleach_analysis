function [fnpt] = DarkBox_modify(image, border, minDist, brightspots, boxwidth, FinalPoints, segmentnumber, outname)
    %% DarkBox Function
    % This function identifies areas within an image that lack fluorescent molecules.
    % It then compiles a list of coordinates (xy) corresponding to the center of the 
    % areas that do not "touch" each other.

    %% Input Parameters:
    % - image: the input image to analyze
    % - border: the border to avoid near the edges of the image
    % - minDist: the minimum distance threshold for spot detection
    % - brightspots: array of detected bright spots
    % - boxwidth: width of the box to search for non-fluorescent areas
    % - FinalPoints: final points array (presumably the result storage)
    % - segmentnumber: a segment identifier
    % - outname: output filename

    % Define output path
    dataout = [pwd '/' outname];
    aIout = image;
    
    % Filter brightspots based on a certain criterion (column 8 > 0)
    alldpts = brightspots(brightspots(:, 8) > 0, :);

    %% Create Grid and Neighborhood
    se = strel('disk', boxwidth, 0);
    sen = se.getnhood();  % Get the neighborhood (structuring element)
    senReal = sen .* 1;
    senReal(boxwidth + 1, boxwidth + 1) = 9;  % Mark the center of the disk
    BGbx = sen - 1;
    BGbx(BGbx == -1) = 1;
    IntNum = sum(sen(:));

    %% Determine Window Size and Create Grid
    window = size(senReal, 1);
    M = ceil(512 / window);
    N = ceil(256 / window);
    newMap = repmat(senReal, M, N);

    %% Find Grid Points
    mrklocalA = [];
    for jk = 1:size(newMap, 1)
        spt = find(newMap(jk, :) == 9)';
        if ~isempty(spt)
            spt(:, 2) = jk;
            mrklocalA = vertcat(mrklocalA, spt);
        end
    end

    %% Filter Grid Points by Border
    posind = mrklocalA;
    posind = posind(posind(:, 1) > border & posind(:, 2) > border & ...
                    posind(:, 2) < 512 - border & posind(:, 1) < 256 - border, :);
    gdbx(:, 1) = posind(:, 1) - boxwidth;
    gdbx(:, 2) = posind(:, 2) - boxwidth;
    gdbx(:, 3) = posind(:, 1) + boxwidth;
    gdbx(:, 4) = posind(:, 2) + boxwidth;
    cnpt = posind;

    %% Calculate Pairwise Distances and Filter by Minimum Distance
    dDistmat = pdist2(cnpt(:, 1:2), alldpts(:, 1:2));  % Pairwise distances between grid points and bright spots
    [blkval, ~] = sort(dDistmat, 2);  % Sort distances in ascending order
    gdbxB = gdbx(blkval(:, 1) > minDist, :);  % Filter based on minimum distance
    cnptB = cnpt(blkval(:, 1) > minDist, :);  % Filter based on minimum distance

    %% Calculate Integrated Intensity for Background Areas
    pts = cnptB;
    tmpBG = [];
    for ff = 1:size(pts, 1)
        pc2 = zeros(1, 6);
        pc2(1, 1) = pts(ff, 1);
        pc2(1, 2) = pts(ff, 2);
        aoix = pts(ff, 2);
        aoiy = pts(ff, 1);

        % Additional processing might go here (the rest is not provided)

    end
    
    %% Visualization
    figure;
    set(gcf, 'position', [696, 29, 637, 952]);
    imagesc(image);
    hold on;


    % Draw boxes around the negative spaces (rejected areas in red, final areas in green)
    for gg = 1:size(bdpt, 1)
        line([bdpt(gg, 1) bdpt(gg, 1)], [bdpt(gg, 2) bdpt(gg, 4)], 'LineWidth', 2, 'Color', 'r');
        line([bdpt(gg, 3) bdpt(gg, 3)], [bdpt(gg, 2) bdpt(gg, 4)], 'LineWidth', 2, 'Color', 'r');
        line([bdpt(gg, 1) bdpt(gg, 3)], [bdpt(gg, 2) bdpt(gg, 2)], 'LineWidth', 2, 'Color', 'r');
        line([bdpt(gg, 1) bdpt(gg, 3)], [bdpt(gg, 4) bdpt(gg, 4)], 'LineWidth', 2, 'Color', 'r');
    end

    for gg = 1:size(fnbx, 1)
        line([fnbx(gg, 1) fnbx(gg, 1)], [fnbx(gg, 2) fnbx(gg, 4)], 'LineWidth', 2, 'Color', 'g');
        line([fnbx(gg, 3) fnbx(gg, 3)], [fnbx(gg, 2) fnbx(gg, 4)], 'LineWidth', 2, 'Color', 'g');
        line([fnbx(gg, 1) fnbx(gg, 3)], [fnbx(gg, 2) fnbx(gg, 2)], 'LineWidth', 2, 'Color', 'g');
        line([fnbx(gg, 1) fnbx(gg, 3)], [fnbx(gg, 4) fnbx(gg, 4)], 'LineWidth', 2, 'Color', 'g');
    end

    %% Finalize and Save the Figure
    movegui center;
    set(gca, 'FontName', 'Arial Black', 'FontSize', 16, 'FontWeight', 'bold', 'LineWidth', 3);
    saveas(gcf, fullfile(dataout, [outname num2str(segmentnumber) '_Negativespace']), 'png');
    hold off;
    close all;
end
