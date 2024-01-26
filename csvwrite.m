function csvwrite(filename, m, r, c, append)
%CSVWRITE Write a comma separated value file.
% Luin Edit: append the data at the end of the file.
%   CSVWRITE(FILENAME,M) writes matrix M into FILENAME as 
%   comma separated values.
%
%   CSVWRITE(FILENAME,M,R,C) writes matrix M starting at offset 
%   row R, and column C in the file.  R and C are zero-based,
%   that is R=C=0 specifies first number in the file.
%
%   See also CSVREAD, DLMREAD, DLMWRITE, WK1READ, WK1WRITE.

%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 5.10.4.3 $  $Date: 2006/11/11 22:44:07 $

%
% test for proper filename
%
if ~ischar(filename)
    error('MATLAB:csvwrite:FileNameMustBeString',...
        'FILENAME must be a string.');
end

%
% Call dlmwrite with a comma as the delimiter
%
if nargin < 3
    r = 0;
end
if nargin < 4
    c = 0;
end
if nargin < 5 || ~strcmpi('-append',append)
    dlmwrite(filename, m, ',', r, c);
else
    dlmwrite(filename, m,'-append', 'delimiter', ',', ...
          'newline', 'pc','roffset', r,'coffset', c);
end;

