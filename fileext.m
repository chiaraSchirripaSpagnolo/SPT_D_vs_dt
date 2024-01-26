function [f,fi]=fileext(type,out2)
% [f,fi]=fileext(type,out2)
% returns a list of file types with their description and eventually a file
% index (if not out2 or for single output) to be used in dialogs concerning
% files (uiputfile, uigetfile, ...)
% fi or f(:,3) contains the file "index", which could be used in variuos
% functions to determine the file type.
% type: 'traj', 'text', 'trackmate', 'bf' (BioFormats importer), 'image',
% 'movie', 'mat'.
% If type ends with '+', an entry for a file mat containing an old file
% list is added, with index=99
% Simplified version without BioFormats importer
if nargin<2||isempty(out2)
    out2=false;
end
%persistent a;
if nargin<1, type='traj';end
if iscell(type)
    f={};
    for typen=1:length(type)
        f=[f;fileext(type{typen})]; %#ok<AGROW>
    end
else
    switch lower(type)
        case {'text','text+'}
            f{1,1} = '*.csv;*.dat;*.txt;*.tsv';f{1,2} = 'Text file (*.csv;*.dat;*.txt;*.tsv)';f{1,3}=1;
            f{2,1} = '*.*';f{2,2} = 'All file (*.*)';f{2,3}=2;
        case {'trackmate','trackmate+'}
            f{2,1} = '*.csv;*.xls';f{2,2} = 'TrackMate spots (*.csv;*.xls)';f{2,3}=12;
            f{1,1} = '*.csv;*.xls';f{1,2} = 'TrackMate spots and tracks (*.csv;*.xls)';f{1,3}=13;
            f{4,1} = '*.*';f{4,2} = 'TrackMate spots (*.*)';f{4,3}=12;
            f{3,1} = '*.*';f{3,2} = 'TrackMate spots and tracks (*.*)';f{3,3}=13;
        case 'traj'
            % try
            %     assert(bfCheckJavaPath(), 'Could not load the Bio-Formats library');
            %     bf=true;
            %     if isempty(a)
            %         a=bfGetFileExtensions();
            %         a=a{1,1};
            %     end
            %     if isempty(a), bf=false; end
            % catch ME
            %    bf=false;
            % end
            f{1,1} = '*.xls';f{1,2} = 'Excel (*.xls)';f{1,3}=1;
            f{2,1} = '*.csv';f{2,2} = 'Comma Separated (*.csv)';f{2,3}=2;
            f{3,1} = '*.dat';f{3,2} = 'Tab Separated (*.dat)';f{3,3}=3;
            f{4,1} = '*.mat';f{4,2} = 'LoadImaris Results (*.mat)';f{4,3}=4;
            f{5,1} = '*.txt';f{5,2} = 'Tab Separated (*.txt)';f{5,3}=5;
            f{6,1} = '*.mat';f{6,2} = '*out.mat -> *TAD*.csv';f{6,3}=6;
            f{7,1} = '*.mat';f{7,2} = 'AnalysisOut.mat -> tracks';f{7,3}=7;
            f{8,1} = '*.mat';f{8,2} = 'u-track output (*.mat)';f{8,3}=8;
            f{9,1} = '*.mat';f{9,2} = 'u-track output, subtrajectories divided (*.mat)';f{9,3}=9;
            % if bf
            %     f{10,1} = a; f{10,2} = 'use u-track (movies)';f{10,3}=10;
            %     f{11,1} = a; f{11,2} = 'use u-track, subtrajectories divided (movies)';f{11,3}=11;
            % end
            f{end+1,1} = '*.mat';f{end,2} = 'old file/dir list';f{end,3}=99;
        case {'image','image+'}
            f={'*.jpg;*.jpeg;*.tif;*.tiff;*.png;*.gif','Image files (*.jpg;*.jpeg;*.tif;*.tiff;*.png;*.gif)',1;...
                '*.bmp;*.pbm;*.pnm;*.pgm;*.ppm','(Bit)map image (*.bmp;*.pbm;*.pnm;*.pgm;*.ppm)',2;...
                '*.jp2;*.jpf;*.jpx;*.j2C;*.j2k','jpeg 2000 (*.jp2;*.jpf;*.jpx;*.j2C;*.j2k)',3; ...
                '*.hdf;*.pcx;*.ras;*.xwd;*.cur;*.ico','Other image format (*.hdf;*.pcx;*.ras;*.xwd;*.cur;*.ico)',4};
        case {'bf','bf+'}
            warning('Luin:fileext','bf not implemented in this simplified version of fileext')
            % try
            %     assert(bfCheckJavaPath(), 'Could not load the Bio-Formats library');
            %     bf=true;
            %     if isempty(a)
            %         a=bfGetFileExtensions();
            %         a=a{1,1};
            %     end
            %     if isempty(a), bf=false; end
            % catch ME
            %     bf=false;
            % end
            % if bf
            %     f{1,1} = a; f{1,2} = 'movies / bioformat import';f{1,3}=10;
            % else
                f=[fileext('image');fileext('movie')];
            % end
        case {'movie','movie+','mov','mov+'}
            f{1,1} = '*.mov;*.mpg;*.avi;*.mj2;*.wmv;*.asf;*.asx;*.mp4;*.m4v';
            f{1,2} = 'movies';
            f{1,3} = 5;
        case {'mat','.mat','mat+','.mat+','*.mat','*.mat+'}
            f{1,1} = '*.mat';
            f{1,2} = 'Old Output (*.mat)';
            f{1,3} = 100;
            
        otherwise
            f{1,1} = '*.*'; f{1,2} = 'All Files (*.*)';f{1,3}=1;
    end
    if isequal(type(end),'+')
        f{end+1,1} = '*.mat';f{end,2} = 'old file/dir list';f{end,3}=99;
%         if size(f,1)>2
%             f{1,3}=11;%???
%         end
    end
end
fi=f(:,3);
if out2 || nargout>1
    f=f(:,1:2);
end
