function para = load_para_yt(para_time)

%ptime = int2str(para_time);

load_str = '/v/raid1a/ytian/MRIdata/Parameters/';

load_str = strcat(load_str,'para_',para_time);

load(load_str)

end