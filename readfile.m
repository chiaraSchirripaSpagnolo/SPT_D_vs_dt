function [output, header,textdata,delm]=readfile(filename,delm,checkall)
%[output, header,textdata,delm]=readfile(filename,delm,checkall)
%if checkall==2, try to use importdata, and return the first line in
%header, numbers in output, textdata like from importdata.
%else if checkall is true, check the length of all the lines. Extremely
%slow.
if nargin>2 && ~isempty(checkall) && checkall==2
    worksheet=importdata(filename);
    header=worksheet.textdata(1,:);
    if isempty(delm),delm=',';end
    if isempty(header{2})
        header=textscan(header{1},'%s','delimiter',delm,'whitespace',' ');
        header=header{1}';
    end
    output=worksheet.data;
    textdata=worksheet.textdata;
%     fid=openw(filename,'rt');
%     d=textscan(fid,'%s','delimiter',delm,'whitespace',' ');
%     fclose(fid);
%     if numel(d{1})==numel(textdata)
        return;
%     end
end
fid=openw(filename,'rt');
str=fgetl(fid);
if nargin<2 || isempty(delm), delm='\t'; end
b=textscan(str,'%s','delimiter',delm,'whitespace',' ');
if length(b{1})==1
    delm=',';
    b=textscan(str,'%s','delimiter',delm,'whitespace',' ');
    if length(b{1})==1
        delm=';';
        b=textscan(str,'%s','delimiter',delm,'whitespace',' ');
        if length(b{1})==1
            delm='\t';
            b=textscan(str,'%s','delimiter',delm,'whitespace',' ');
        end
    end
end
if str(end)==delm
    b{1}(end+1)={''};
end
% b=textscan(str,'%s','delimiter','\t','whitespace',' ','BufSize',65535);
lngth=size(b{1},1);
if nargin<3 || isempty(checkall) || ~checkall
    d=textscan(fid,'%s','delimiter',delm,'whitespace',' ');
    % d=textscan(fid,'%s','delimiter','\t','whitespace',' ','BufSize',65535);
    e=reshape(d{1}(1:lngth*floor(size(d{1},1)/lngth)),lngth,floor(size(d{1},1)/lngth))';
    checkall=0;
else
    e=b{1}';
    while(~isequal(str,-1))
        str=fgetl(fid);
        if ~isequal(str,-1)
            if isempty(str)
                b={{''}};
            else
                b=textscan(str,'%s','delimiter',delm,'whitespace',' ');
                if ~iscell(b{1}), b={b}; end
            end
            diffl=length(b{1})-lngth;
            if diffl==0
                e=[e;b{1}'];
            elseif diffl<0
                e=[e;b{1}',repmat({''},1,-diffl)];
            else
                lngth=length(b{1});
                e=[e,repmat({''},size(e,1),diffl);b{1}']; %#ok<*AGROW>
            end
        end
    end
    b={e(1,:)'};
    e(1,:)=[];
end
%floor to take care of empty line(s) at the end of file.
switch nargout
    case 1
        output=[b{1}';e];
    case 2
        output=str2double(e);
        header=b{1}';
    otherwise
%         output=str2double(e);
%         frstlout=str2double(b{1}');
%         if ~all(isnan(frslout))
%             output=[frstlout;e];
%         end
        if ~checkall
            output=str2double(e);
            frstlout=str2double(b{1}');
            if ~all(isnan(frstlout))
                output=[frstlout;e];
            end
        else
            output=str2double([b{1}';e]);
        end
        header=b{1}';
        textdata=[b{1}';e];
end
fclose(fid);
end