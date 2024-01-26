optdlg.Resize='on';
optdlg.WindowStyle='normal';
optdlg.Interpreter='none';
%     root='prova';
%     count=0;
%     fpath=pwd;
%     nsubmin=+Inf;
%     div=true;
%     nmin=0;
%     tmin=0;
%     kadd=0;
%     tsigma=2.58;
%     onlyending=true;
%     ampRatioLimit=1;
if pref('is','LuosTrack','dividepref')
    answ=pref('get','LuosTrack','dividepref');
else
    defs={'false', 'true','+Inf','2.58','1','true'};
    prompt={'If applicable, correct merge and split? (0/1)',...
        'Consider only terminal or initial subtrajectories? (0/1)',...
        'Max number of spots for considering unlinking of subtrajectory',...
        'Consider uncertainties with coefficients (n of sigma):',...
        'Max intensities ratio for merge and split',...
        'Apply this preferences to all trajectories?'} ;
%         'Max intensities ratio for joining',...
%         
    answ=inputdlg(prompt,'Separate wrong "merged/splitted" subtrajectories',1,defs,optdlg);
    drawnow; pause(0.05);
    if ~isempty(str2num(answ{end})) && str2num(answ{end})
        pref('add','LuosTrack','dividepref',answ(1:end-1));
    end
    answ(end)=[];
end
div=str2num(answ{1}); %#ok<*ST2NM>
onlyending=str2num(answ{2});
nsubmin=str2num(answ{3});
tsigma=str2num(answ{4});
ampRatioLimit=str2num(answ{5});
if isempty(div),div=false;end
if isempty(onlyending),onlyending=false;end

