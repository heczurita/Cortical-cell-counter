function calculateCellsPositionsMultipleImagesGraph()
    % Select images
    [fileNames, pathName] = uigetfile(...
        {'*.jpg;*.png;*.jpeg;*.tif;*.tiff','Image Files (JPEG, PNG, TIFF)'},...
        'Select Images of the Cortical Columns','MultiSelect','on');
    if isequal(fileNames,0); disp('User canceled'); return; end
    if ischar(fileNames); fileNames = {fileNames}; end

    % Select channels
    channelOptions = {'Green','Red','Blue'};
    [selCh, okCh] = listdlg(...
        'PromptString','Select channels present:',...
        'SelectionMode','multiple','ListString',channelOptions);
    if ~okCh; disp('No channels selected'); return; end
    channelLabels = channelOptions(selCh);
    channelFields = cellfun(@matlab.lang.makeValidName, channelLabels, 'UniformOutput', false);

    % Ask if user wants to delineate layers
    layerDesc = questdlg('Delineate custom layers?','Layers','Yes','No','Yes');
    if strcmp(layerDesc,'Yes')
        layerOptions = {'L1','L2','L3','L2/3','L4','L5','L5a','L5b','L6a','L6b','L6'};
        [selLy, okLy] = listdlg(...
            'PromptString','Select layers to delineate:',...
            'SelectionMode','multiple','ListString',layerOptions);
        if okLy && ~isempty(selLy)
            selectedLayers = layerOptions(selLy);
            nLayers = numel(selectedLayers);
            boundaries = zeros(nLayers,2);
            for j = 1:nLayers
                prompt = {sprintf('Enter top and bottom normalized boundaries for %s (format: top,bottom)', selectedLayers{j})};
                def = {'1,0'};
                answer = inputdlg(prompt, sprintf('Boundary for %s', selectedLayers{j}), [1 50], def);
                if isempty(answer); disp('Layer input canceled'); return; end
                nums = sscanf(answer{1}, '%f,%f'); boundaries(j,:) = nums';
            end
            [~, idxSort] = sort(boundaries(:,1), 'descend');
            edgesDesc = [boundaries(idxSort,1); boundaries(idxSort(end),2)];
            layerNamesDesc = selectedLayers(idxSort);
        else
            edgesDesc = [];
            layerNamesDesc = {};
        end
    else
        edgesDesc = [];
        layerNamesDesc = {};
    end

    % Initialize storage
    cellPositions = struct();
    summaryData = cell(1, numel(fileNames));

    % Define dark/light colors for pie charts
    dark = struct('Green', [0,0.5,0], 'Red', [0.5,0,0], 'Blue', [0,0,0.5]);
    light = struct('Green', [1,1,0], 'Red', [1,0.5,0], 'Blue', [0,1,1]);

    % Process each image
    for i = 1:numel(fileNames)
        safeName = matlab.lang.makeValidName(fileNames{i});
        img = imread(fullfile(pathName, fileNames{i}));
        if size(img,3) == 1; img = repmat(img, [1,1,3]); end

        % Select column bounds
        h = figure('Name', fileNames{i}); imshow(img);
        title('Select Top and Bottom of Column');
        [~, y] = ginput(2); y = sort(y); close(h);
        yTop = y(1); yBottom = y(2);

        ratiosByChannel = struct();
        % Per-channel click collection with inverted contrast
        for c = 1:numel(channelLabels)
            lbl = channelLabels{c}; fld = channelFields{c};
            chImg = selectChannelImage(img, lbl);
            % Invert image
            if isinteger(chImg)
                invImg = imcomplement(chImg);
            else
                invImg = 1 - mat2gray(chImg);
            end

            hCh = figure('Name', sprintf('%s - %s (inverted)', fileNames{i}, lbl));
            imshow(invImg); colormap(gca, 'gray'); hold on;
            title(sprintf('Left-click add, right-click undo (%s). Other to finish', lbl));
            plot([1 size(img,2)], [yTop yTop], 'g-', 'LineWidth', 1);
            plot([1 size(img,2)], [yBottom yBottom], 'g-', 'LineWidth', 1);

            xPts = []; yPts = [];
            while true
                [x, y, btn] = ginput(1);
                if isempty(btn) || ~ismember(btn, [1,3]); break; end
                if btn == 1 % add
                    xPts(end+1) = x; yPts(end+1) = y;
                    plot(x, y, 'o', 'MarkerEdgeColor', lbl, 'MarkerFaceColor', lbl, 'MarkerSize', 6);
                else % undo
                    if ~isempty(xPts)
                        xPts(end) = []; yPts(end) = [];
                        cla; imshow(invImg); hold on;
                        plot([1 size(img,2)], [yTop yTop], 'g-', 'LineWidth', 1);
                        plot([1 size(img,2)], [yBottom yBottom], 'g-', 'LineWidth', 1);
                        for k = 1:numel(xPts)
                            plot(xPts(k), yPts(k), 'o', 'MarkerEdgeColor', lbl, 'MarkerFaceColor', lbl, 'MarkerSize', 6);
                        end
                    end
                end
            end
            close(hCh);

            ratiosByChannel.(fld) = (yPts - yTop) / (yBottom - yTop);
            cellPositions.(safeName).(fld) = ratiosByChannel.(fld);
        end
        summaryData{i} = struct('Image', fileNames{i}, 'Labels', {channelLabels}, 'Fields', {channelFields}, 'Ratios', ratiosByChannel);
    end

    % Figure 1: Scatter
    plotScatterWithLayers(summaryData, channelLabels, channelFields, edgesDesc, layerNamesDesc);
    % Figure 2: Pie-by-layer if layers defined
    if ~isempty(edgesDesc)
        plotPieByLayer(summaryData, channelLabels, channelFields, edgesDesc, layerNamesDesc, dark, light);
    end
    % Figure 3: Distribution of each channel
    plotDistributionByChannel(summaryData, channelLabels, channelFields);

    % Save results to workspace and Excel
    assignin('base', 'cellPositions', cellPositions);
    disp('cellPositions (normalized ratios) saved to workspace');

    excelFile = fullfile(pathName, 'CellPositions.xlsx');
    T = table();
    for i = 1:numel(fileNames)
        fname = fileNames{i}; safeName = matlab.lang.makeValidName(fname);
        for f = fieldnames(cellPositions.(safeName))'
            ratios = cellPositions.(safeName).(f{1});
            for idx = 1:numel(ratios)
                T = [T; {fname, f{1}, idx, ratios(idx)}];
            end
        end
    end
    T.Properties.VariableNames = {'Image','Channel','Index','Ratio'};
    writetable(T, excelFile);
    disp(['Excel saved to ', excelFile]);
end

function chImg = selectChannelImage(img, ch)
    switch ch
        case 'Green', chImg = img(:,:,2);
        case 'Red',   chImg = img(:,:,1);
        case 'Blue',  chImg = img(:,:,3);
    end
end

function plotScatterWithLayers(data, labels, fields, edgesDesc, layerNamesDesc)
    colors = {'g','r','b'}; numCh = numel(fields); offsets = 0:(numCh-1);
    figure('Name','Figure 1 - Cell Positions by Channel'); hold on;
    if ~isempty(edgesDesc)
        for k = 2:numel(edgesDesc)-1, yline(edgesDesc(k), 'r:', 'LineWidth', 2); end
    end
    for c = 1:numCh
        arr = [];
        for idx = 1:numel(data)
            arr = [arr; data{idx}.Ratios.(fields{c})(:)];
        end
        xJ = offsets(c) + (rand(size(arr))*0.6 - 0.3);
        scatter(xJ, arr, 36, colors{c}, 'filled','MarkerFaceAlpha', 0.6);
        text(offsets(c), 0.97, sprintf('n=%d', numel(arr)), 'Color', colors{c}, 'HorizontalAlignment', 'center');
    end
    xlim([-0.5, offsets(end)+0.5]); xticks(offsets); xticklabels(labels);
    ylim([0 1]); ylabel('Normalized Position'); xlabel('Channel');
    if ~isempty(edgesDesc)
        yyaxis right;
        mids = (edgesDesc(1:end-1) + edgesDesc(2:end)) / 2;
        [ms, idx] = sort(mids);
        yticks(ms); yticklabels(layerNamesDesc(idx)); ylabel('Layer'); set(gca, 'YColor', 'k'); yyaxis left;
    end
    hold off;
end

function plotPieByLayer(data, labels, fields, edgesDesc, layerNamesDesc, dark, light)
    figure('Name','Figure 2 - Layer Pie Charts');
    edgesAsc = sort(edgesDesc);
    layerNamesAsc = flipud(layerNamesDesc(:));
    for c = 1:numel(fields)
        arr = [];
        for idx = 1:numel(data)
            arr = [arr; data{idx}.Ratios.(fields{c})(:)];
        end
        counts = histcounts(arr, edgesAsc);
        props = counts / sum(counts);
        nz = props > 0; propsNZ = props(nz); namesNZ = layerNamesAsc(nz);
        dColor = dark.(labels{c}); lColor = light.(labels{c});
        m = numel(propsNZ); t = linspace(0,1,m)'; cmap = (1-t)*dColor + t*lColor;
        subplot(1, numel(fields), c);
        h = pie(propsNZ); patches = findobj(h, 'Type', 'patch');
        for j = 1:numel(patches), set(patches(j),'FaceColor', cmap(j,:)); end
        title(sprintf('%s Channel', labels{c})); legend(namesNZ,'Location','eastoutside');
    end
end

function plotDistributionByChannel(data, labels, fields)
    figure('Name','Figure 3 - Distribution of Channels'); hold on;
    % Specify colors for each channel
    colorMap = struct('Green', 'g', 'Red', 'r', 'Blue', 'b');
    for c = 1:numel(fields)
        arr = [];
        for idx = 1:numel(data)
            arr = [arr; data{idx}.Ratios.(fields{c})(:)];
        end
        [f, xi] = ksdensity(arr);
        chColor = colorMap.(labels{c});
        plot(xi, f, 'Color', chColor, 'LineWidth', 1.5);
    end
    legend(labels, 'Location', 'best');
    xlabel('Normalized Position'); ylabel('Density');
    title('Channel Distribution');
    xlim([0 1])
    hold off;
end
