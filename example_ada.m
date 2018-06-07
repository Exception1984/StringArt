clc;
clear;
close all;

% Ada Lovelace

% INFO: Pass an importance map using 'importanceMapPath', '<path>', ...
stringArt( ...
    'inputFilePath', 'input/ada.png', ... % path to the input file
    'outputDirPath', 'output/ada', ... % path to the output directory
    'outputFileNamePrefix', 'ada', ... % prefix of the output files
    'invertInput', true, ... % false -> reconstruct white area, true -> reconstruct black area
    'invertOutput', true, ... % false -> white string, true -> black string
    'numPins', 256, ... % number of pins that are placed on the frame (default is 256)
    'threadThickness', 0.15, ... % physical thickness of the thread in mm (default is 0.15)
    'frameDiameter', 614.4, ... % physical diameter of the circular frame in mm (default is 614.4)
    'pinSideLength', 2, ... % physical side length of a pin with quadratic cross section in mm (default is 2)
    'superSamplingWindowWidth', 8, ... % side length of the super sampling window (default is 8)
    'minAngle', pi/8, ...
    'dataPath', './data'); % Minimum angle (measured from frame center) between two connected pins (default is pi/8)