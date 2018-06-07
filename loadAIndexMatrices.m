function [EPI, EPV, CEI, ITI] = loadAIndexMatrices(matrixPath, imSize, invalidEdgeIndices)  
    CEI = matrixFromFile([matrixPath '/A_s=' num2str(imSize) '_CEI'], 'uint32');
    edgeMask = true(size(CEI));
    
    if nargin > 3
        edgeMask = ~ismember(CEI, invalidEdgeIndices);
    end
    
    CEI = CEI(edgeMask);
    
    EPI = matrixFromFile([matrixPath '/A_s=' num2str(imSize) '_EPI'], 'uint32');
    EPI = EPI(edgeMask);
    
    EPV = matrixFromFile([matrixPath '/A_s=' num2str(imSize) '_EPV'], 'double');
    EPV = EPV(edgeMask);

    if nargout > 3
        % Load ITI Matrix
        % Read parameters
        fileID = fopen([matrixPath '/A_s=' num2str(imSize) '_ITI.txt']);
        params = textscan(fileID, '%s %s %d');
        params = params{3};
        n = double(params(2));
        fclose(fileID);

        m = numel(EPI);

        ITI = sparse((1 : m)', double(EPI), true(m, 1), m, n);
    end
end