function [cardiac_signal,resp_signal,mask] = get_cardiac_signal_SPGR(Image)
image_for_std = crop_half_FOV(abs(Image(:,:,1:28,:)));
nslice = size(image_for_std,4);
sx = size(image_for_std,1);
nof = size(Image,3);
image_for_std = reshape(image_for_std,[sx,sx,7,4,nslice]);
std_map = std(image_for_std,0,3);
std_map = squeeze(sum(std_map,4));
std_map = imgaussfilt(std_map,5);
filter = fspecial('gaussian',sx,sx/10);
std_map = std_map.*filter;
mask = zeros(sx,sx,nslice);
for i=1:nslice
    [x,y] = find(std_map(:,:,i)==max(max(std_map(:,:,i))));
    mask(x,y,i) = 1;
    mask(:,:,i) = bwdist(mask(:,:,i));
end

mask = mask<45;
% mask_resp = sum(mask,1);
% mask_resp = mask_resp==0;



dt = 35;
Fs = 1000/dt;
df = Fs/nof;


% cardiac_range = round((0.5/df):(2.2/df));

fpass_cardiac = [0.5/df/(nof/2),2.2/df/(nof/2)];


cardiac_signal = permute(mask,[1,2,4,3]).*crop_half_FOV(abs(Image(:,:,:,:)));
cardiac_signal = permute(cardiac_signal,[1,2,4,3]);
cardiac_signal = reshape(cardiac_signal,[sx*sx*nslice,nof]);
cardiac_signal = sum(cardiac_signal,1);
cardiac_signal = bandpass(cardiac_signal,fpass_cardiac);
% cardiac_signal(cardiac_signal(:,1)==0,:) = [];
% cardiac_coeff = pca(cardiac_signal);
% cardiac_signal_fft = fft(cardiac_coeff(:,1:10)./sos(cardiac_coeff(:,1:10),1),[],1);
% signal_fft = fftshift(signal_fft,1);

% [~,idx] = max(max(abs(cardiac_signal_fft(cardiac_range,:))),[],2);
% [~,idx_max] = max(abs(cardiac_signal_fft(cardiac_range,idx)));
% idx_max = idx_max + cardiac_range(1)-1;

% cardiac_signal = cardiac_coeff(:,idx);

% mask_resp = false(size(mask));
% 
% for i=1:nslice
%     [x,y] = find(mask(:,:,i));
%     left(i) = min(y);
%     right(i) = max(y);
%     up(i) = min(x);
%     down(i) = max(x);
%     mask_resp(up:down,1:left,i) = true;
% end

resp_range = round((0.2/df):(0.5/df));
fpass_resp = [0.2/df/(nof/2),0.5/df/(nof/2)];

resp_signal = squeeze(sum(abs(Image(:,:,:,:)),2));
resp_signal = permute(resp_signal,[1,3,2]);
resp_signal = reshape(resp_signal,[sx*2*nslice,nof]);
resp_coeff = pca(resp_signal);
resp_signal_fft = fft(resp_coeff(:,1:10),[],1);

[~,idx] = max(max(abs(resp_signal_fft(resp_range,:))),[],2);

[~,idx_max] = max(abs(resp_signal_fft(resp_range,idx)),[],1);
peak_frequency = (resp_range(1) + idx_max - 1)*df/Fs;

resp_signal = lowpass(resp_coeff(:,idx),peak_frequency);

% 
% resp_signal = permute(mask_resp,[1,2,4,3]).*crop_half_FOV(abs(Image(:,:,:,:)));
% resp_signal = permute(resp_signal,[1,2,4,3]);
% resp_signal = reshape(resp_signal,[sx*sx*nslice,nof]);
% resp_signal = sum(resp_signal,1);
% resp_signal = bandpass(resp_signal,fpass_resp);
% 
% resp_signal(resp_signal(:,1)==0,:) = [];
% resp_coeff = pca(resp_signal);
% resp_signal_fft = fft(resp_coeff(:,1:5)./sos(resp_coeff(:,1:5),1),[],1);
% 
% 
% [~,idx] = max(max(abs(resp_signal_fft(resp_range,:))),[],2);
% resp_signal = bandpass(resp_coeff(:,idx),fpass_resp);

% cardiac_signal = sum(sum(sum(cardiac_signal)),4);
% cardiac_signal = squeeze(cardiac_signal);
% 
% resp_signal = crop_half_FOV(abs(Image(:,:,:,:))).*permute(mask_resp,[1,2,4,3]);

end