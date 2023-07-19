function [save_dir_time] = save_para_yt(para)

para_temp = para;
save_dir = '/v/raid1a/ytian/MRIdata/Parameters/';
para_name = dir(save_dir);
load(strcat(save_dir,para_name(end).name));
try
    setn = para.set;
end
clear para
para = para_temp;
if exist('setn','var')
    para.set = setn+1;
else
    para.set = 0;
end
save_dir_time =  datestr(clock,'yymmdd_hhMMSS');
para.time = save_dir_time;
MID = strfind(para.dir.load_kSpace_name,'MID');
MID = para.dir.load_kSpace_name(MID+5:MID+7);
para.MID = MID;
para.result_name = strcat('MID_',MID,'_',save_dir_time);
save_name = strcat('para_',save_dir_time);

save([save_dir,save_name],'para')
disp(para.time)
disp(strcat('set',num2str(setn+1)))

return