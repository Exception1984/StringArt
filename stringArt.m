function [A, x] = stringArt(varargin)
    % Set Default Values
    invertInput = true;
    invertOutput = true;
    numPins = 256;
    threadThickness = 0.15;
    frameDiameter = 614.4;
    pinSideLength = 2;
    superSamplingWindowWidth = 8;
    minAngle = 0;
    importanceMapPath = '';
    range = 1.0;
    dataPath = './data';
    
    % Parse Input Parameters
    for i = 1 : 2 : numel(varargin)
        switch varargin{i}
            case 'inputFilePath'
                fileName = varargin{i + 1};
            case 'outputDirPath'
                outputPath = varargin{i + 1};
            case 'outputFileNamePrefix'
                outputFileNamePrefix = varargin{i + 1};
            case 'invertInput'
                invertInput = varargin{i + 1};
            case 'invertOutput'
                invertOutput = varargin{i + 1};
            case 'numPins'
                numPins = varargin{i + 1};
            case 'threadThickness'
                threadThickness = varargin{i + 1};
            case 'frameDiameter'
                frameDiameter = varargin{i + 1};
            case 'pinSideLength'
                pinSideLength = varargin{i + 1};
            case 'superSamplingWindowWidth'
                superSamplingWindowWidth = varargin{i + 1};
            case 'minAngle'
                minAngle = varargin{i + 1};
            case 'importanceMapPath'
                importanceMapPath = varargin{i + 1};
            case 'dataPath'
                dataPath = varargin{i + 1};
        end
    end
    
    % Super Sampling Window Width must be power of 2
    assert(mod(superSamplingWindowWidth, 2) == 0);
    
    % Adapt Thread Thickness to allow a domain width that is power of 2
    threadThickness = adaptThreadThickness(threadThickness, frameDiameter);

    imOrig = double(imread(fileName)) / 255;
    
    highRes = round(frameDiameter / threadThickness);
    lowRes = highRes / superSamplingWindowWidth;

    %Resize image
    if size(imOrig, 1) ~= lowRes
        imOrig = imresize(imOrig, [lowRes, lowRes]);
    end
    
    % Stretch Histogram of input image to enhance contrast
    imOrig = stretchHistogram(imOrig);
    
    % Shorten range of input image to remove pure white regions
    imOrig = imOrig * range;
    
    img = imOrig;
    
	targetFilePath = [outputPath '/' outputFileNamePrefix '-target.png'];
    
    [targetImage, mask] = maskImage(imOrig, 0.5 * size(imOrig, 1));
    
	if xor(invertInput, invertOutput) == true
        targetImage = 1 - targetImage;
	end
    
    targetImage(mask) = 1;
    
    % Write Target Image
    filepath = fileparts(targetFilePath);
    if ~(exist(filepath, 'dir') == 7)
        mkdir(filepath);
    end
    imwrite(targetImage, targetFilePath);
    
    mask = ~mask;
    
    if (invertInput)
        img = imcomplement(img);
    end
    
	importanceMap = ones(size(img));
    
    % Flip image upside down to accout for different coordinate systems of
    % Matlab (left-handed) and robot (right-handed)
    img = flipud(img);
    
    % Load Importance Map if provided
    if ~strcmp(importanceMapPath, '')
        importanceMap = double(imread(importanceMapPath)) / 255;
        
        if size(importanceMap, 1) ~= size(img, 1) || ...
            size(importanceMap, 2) ~= size(img, 2)
            importanceMap = imresize(importanceMap, size(img));
        end
        
        importanceMap = flipud(importanceMap);
    end
    
    % Precompute Matrix
    matrixPath = precomputeMatrix(numPins, threadThickness, frameDiameter, pinSideLength, superSamplingWindowWidth, dataPath);
    
	pickedEdgesSequencePath = [outputPath '/' outputFileNamePrefix '_NH-' num2str(numPins) '_DW-' num2str(lowRes) '_WS-' num2str(superSamplingWindowWidth) '.txt' ];
    
    % Check for existence of file
    if exist(pickedEdgesSequencePath, 'file') == 2
        e = edgeCodes(numPins);
        x = zeros(size(e, 1), 1);
        
        fileID = fopen(pickedEdgesSequencePath, 'r');
        pickedEdgesSequence = fscanf(fileID, '%d');
        fclose(fileID);
        
        x(pickedEdgesSequence) = 1;
    else
        [x, pickedEdgesSequence] = optimizeStringsGreedyMultiSampling(img, lowRes, superSamplingWindowWidth, minAngle, numPins, importanceMap, matrixPath);
        
        % Write EdgeList to file
        fileID = fopen(pickedEdgesSequencePath, 'w');
        fprintf(fileID, '%i\n', pickedEdgesSequence);
        fclose(fileID);

        %Restore imOrig
        lowResImageWidth  = lowRes;
        highResImageWidth = lowRes * superSamplingWindowWidth;

        A = loadFullMatrix(matrixPath, highResImageWidth);

        recon = min(1, A * x);
        reconImageHigh = reshape(recon, [highResImageWidth, highResImageWidth]);
        reconImageHigh = flipud(reconImageHigh');

        reconImageLow = imresize(reconImageHigh, [lowResImageWidth, lowResImageWidth]);

        if (invertOutput)
            reconImageHigh = 1 - reconImageHigh;
            reconImageLow = 1 - reconImageLow;
        end

        recon512Path = [outputPath '/' outputFileNamePrefix '_NH-' num2str(numPins) '_DW-' num2str(lowRes) '_WS-' num2str(superSamplingWindowWidth) '-LowRes.png'];
        reconNativePath = [outputPath '/' outputFileNamePrefix '_NH-' num2str(numPins) '_DW-' num2str(lowRes) '_WS-' num2str(superSamplingWindowWidth) '-HighRes.png'];

        imwrite(reconImageHigh, reconNativePath);
        imwrite(reconImageLow, recon512Path);

        rmseValue = rmse(reconImageLow(mask), targetImage(mask));

        % Write RMSE File
        fileID = fopen([outputPath '/' outputFileNamePrefix '_RMSE-Value.txt' ], 'w');
        fprintf(fileID, '%16.16f\n', rmseValue);
        fclose(fileID);

        % Write Error Image
        gradient = ColorGradient(0.0, 1.0);
        errorImage = abs(targetImage - reconImageLow);
        heatMap = gradient.ColorAtValue(errorImage);
        imwrite(heatMap, [outputPath '/' outputFileNamePrefix '_ErrorImage.png']);
    end
    
    % Compute a consecutive Path and write Result File
    resultFilePath = [outputPath '/' outputFileNamePrefix '_NH-' num2str(numPins) '_DW-' num2str(lowRes) '_WS-' num2str(superSamplingWindowWidth) '_RESULT.txt' ];
    outputToConsecutive(x, resultFilePath, matrixPath);
end