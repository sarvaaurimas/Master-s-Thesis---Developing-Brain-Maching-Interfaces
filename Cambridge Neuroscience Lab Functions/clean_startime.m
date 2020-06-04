function start_time = clean_startime(sync_mess)
fileID = fopen(sync_mess);
C = textscan(fileID,'%s');
ss=(C{1}(13));
tim_tmp = cell2mat(ss);
start_time = str2num(tim_tmp(1:(strfind(tim_tmp,'@')-1)));
end