function T1_map = SASHA_fitting(im,kSpace_info)
nor = 120;



siz = size(im);
if length(siz)<4
    siz(4) = 1;
end

im = abs(im);
PD = im(:,:,1,:);
im = im(:,:,2:end,:);
PD(PD<1) = 1;
im_normalized = im./PD;

SRTs(1) = kSpace_info.Protocol.alTI/1000;
for i=1:siz(3)-2
    SRTs(end+1) = SRTs(end) + kSpace_info.Protocol.sWipMemBlock.alFree(4);
end



T1 = 1:3500;
flip_angle = kSpace_info.Protocol.adFlipAngleDegree;
PD_flip_angle = 2;
TR = kSpace_info.Protocol.alTR(2)/1000;
nor = kSpace_info.Protocol.sKSpace.lRadialViews;

dic = Bloch_SASHA(SRTs,T1,flip_angle,nor,TR,PD_flip_angle);


% dic_norm = sos(dic,2);
% dic = dic./dic_norm;

im_normalized = permute(im_normalized,[5,3,1,2,4]);
im_normalized = reshape(im_normalized,[1,siz(3)-1,siz(1)*siz(2)*siz(4)]);
% im_norm = sos(im_normalized,2);
% im_normalized = im_normalized./im_norm;

dic = mean(dic,2);
im_normalized = mean(im_normalized,2);

d = im_normalized-dic;
d = sos(d,2);
[~,T1_map] = min(d,[],1);
T1_map = reshape(T1_map,[siz(1),siz(2),siz(4)]);

