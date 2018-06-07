function writeResultFile(stringList, outsideMask, filePath)
    fileID = fopen(filePath, 'w');
    
    data = [stringList double(outsideMask)];
    
    fprintf(fileID, 'FROM-Pin\tTO-Pin\tOUT-Edgetype\tIN-Edgetype\tOUTSIDE\n');
    fprintf(fileID, '%d\t%d\t%d\t%d\t%d\n', data');
    fclose(fileID);
end