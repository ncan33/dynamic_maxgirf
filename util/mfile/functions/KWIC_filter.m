function kwic = KWIC_filter(sx,nor,center)

kwic = zeros(sx,nor,'single');
kwic(:,1) = 1;
tan_theta = sx/2/nor;
for i=1:nor
    kwic(1:round((nor-i)*tan_theta),i) = 1;
end
kwic(sx/2+1:end,:) = flipud(kwic(1:sx/2,:));
kwic(:,2:3:end) = kwic(:,1:3:end);
kwic(:,3:3:end) = kwic(:,1:3:end);

if center<nor/2
    kwic = circshift(kwic,[0,center]);
    kwic(:,1:center) = fliplr(kwic(:,(1:center)+3+center));
else
    center = nor-center;
    kwic = circshift(kwic,[0,center]);
    kwic(:,1:center) = fliplr(kwic(:,(1:center)+3+center));
    kwic = fliplr(kwic);
end
%figure,imagesc(kwic)