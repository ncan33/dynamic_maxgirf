function filter_all = ramp_filter_for_pre_interp_3D(para)
nos = length(para.Recon.nor);% number of slice
[X,Y] = meshgrid(1:para.Recon.sx,1:para.Recon.sx);
X = X - (para.Recon.sx+1)/2;
Y = Y - (para.Recon.sx+1)/2;
filter = sqrt(X.^2 + Y.^2);
fully_sampled_radius = para.Recon.nor*2/pi;
filter_all = repmat(filter,[1 1 nos]);
for i=1:nos
    filter_temp = filter;
    filter_temp(filter_temp<fully_sampled_radius(i)) = fully_sampled_radius(i);
    filter_all(:,:,i) = filter_temp;
end

filter_all = filter_all/(para.Recon.sx+1)*2;
filter_all(filter_all>1) = 1;
filter_all = fftshift3(filter_all);
filter_all = single(filter_all);

end