function T1 = double_FA_T1(ref_image,FA,TR)

% FA = FA/180*pi;
% FA_PD = 2/180*pi;

siz = size(ref_image);
% ref_image = ref_image(:);
% ref_image = ref_image.*mask;

% B1(B1<0.5) = 0.5;
% B1(B1>1.3) = 1.3;
% B1 = B1.*mask;

ref_image = reshape(ref_image,[siz(1)*siz(2),siz(3)]);
% B1 = reshape(B1,[siz(1)*siz(2),siz(3)]);

% T1 = zeros([prod(siz),1]);

% num = ref_image.*sin(B1*FA_PD) - sin(B1*FA);
% den = ref_image.*cos(B1*FA).*sin(B1*FA_PD) - cos(B1*FA_PD).*sin(B1*FA);
% E1 = num./den;
% 
% T1 = 1./abs(-log(E1)/TR);
% T1 = reshape(T1,siz);
% T1(T1>3000) = 3000;

% B1 = round(B1*100)-49;

T1 = zeros([siz(1)*siz(2),siz(3)]);

dic_all = Quant.get_dic_all(FA,TR,1);

for nslice=1:siz(3)
    tic
    im_temp = ref_image(:,nslice);
    dic_temp = dic_all(nslice,:,:);
    d = im_temp - dic_temp;
    d = abs(d);
    [~,T1_temp] = min(d,[],2);
    T1(:,nslice) = T1_temp;
%     B1_temp = B1(:,nslice);
%     for i=1:siz(1)*siz(2)
%         dic_temp_temp = dic_temp(:,:,B1_temp(i));
%         d = im_temp(i) - dic_temp_temp;
%         d = abs(d);
%         [~,T1_temp] = min(d,[],2);
%         T1(i,nslice) = T1_temp;
%     end
    toc
end

% 
% for i=1:prod(siz(1:2))
%     fprintf([num2str(i/prod(siz(1:2))),'\n'])
%     if B1(i)~=0
%         dic = Quant.get_dictionary_with_phase_interleaving_3_SMS_3(FA,TR,B1(i));
%         for j=1:siz(3)
%             d = abs(ref_image(i,j) - dic(j,:));
%             [~,T1(i,j)] = min(d,[],2);
%         end
%     end
% end


T1 = reshape(T1,siz);