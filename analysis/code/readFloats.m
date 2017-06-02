function [dat] = readFloats(filename,ncol)
%Read a file, filename, with ncol columns of floats into a cell array. If it's a csv,
%designate the comma delimiter.
    
    [~,~,extension] = fileparts(filename);
    floatString=repmat('%f',1,ncol);
    
    fid = fopen(filename);
    if strcmp(extension,'.csv')
        dat=cell2mat(textscan(fid,floatString,'Delimiter',','));
    else
        dat=cell2mat(textscan(fid,floatString));
    end
    fclose(fid);

end

