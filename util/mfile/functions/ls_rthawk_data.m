function ls_rthawk_data(area, vol_number, file_index)
if nargin == 2
    file_index = [];
end

if ismac
    path_server = dir(sprintf('/Users/ytian/MRIdata/%s/*%d*', area, vol_number));
    path_server = fullfile(path_server.folder, path_server.name);
    path_yt = path_server;
else
    if contains(area, 'phantom')
        path_server = dir(sprintf('/server/sdata/ytian/mri_data/disc/phantom/*%d*', vol_number));
    else
        path_server = dir(sprintf('/server/sdata/mri_data/disc/%s/vol%04g_*', area, vol_number));
    end
    path_server = fullfile(path_server.folder, path_server.name);
    
    path_yt = strfind(path_server, 'mri_data');
    path_yt = path_server(path_yt:end);
    path_yt = ['/server/sdata/ytian/', path_yt];
end


if ~isfolder(path_yt)
    mkdir(path_yt)
end

path_yt_raw = [path_yt, '/raw_hawk/'];

if ~isfolder(path_yt_raw)
    mkdir(path_yt_raw)
end

path_server_raw = [path_server, '/raw_hawk/'];

path_search_all = {};
path_search_all{1} = path_server_raw;
sub_folder_temp = search_sub_folder(path_search_all{1});
sub_folder_to_be_searched = sub_folder_temp;
sub_folder_all = [path_server_raw; sub_folder_temp];
while ~isempty(sub_folder_to_be_searched)
    sub_folder_temp = search_sub_folder(sub_folder_to_be_searched{1});
    sub_folder_to_be_searched(1) = [];
    sub_folder_to_be_searched = [sub_folder_to_be_searched; sub_folder_temp];
    sub_folder_all = [sub_folder_all; sub_folder_temp];
end

all_dat = [];
for K = 1 : length(sub_folder_all)
    thisdir = sub_folder_all{K};
    dat_temp = dir(fullfile(thisdir, 'usc_disc_yt*raw.dat'));
    all_dat = [all_dat; dat_temp];
end

% all_dat = dir([path_server_raw, '*raw.dat']);

n_files = length(all_dat);
if isempty(file_index)
    file_index = 1:n_files;
end

for i = file_index %1:n_files
        file_name = fullfile(all_dat(i).folder, all_dat(i).name);
        fprintf('--------------------------------- \n')
        fprintf('data number:  %d/%d \n', i, n_files)
        fprintf('file name:     %s \n',file_name)
end
end