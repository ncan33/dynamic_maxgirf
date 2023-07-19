function Image_patch = get_patch(Image,patch_size)

[sx,sy,nof] = size(Image);
px = patch_size(1);
py = patch_size(2);

Nx = sx/px;
Ny = sy/py;

% Image_patch = zeros(px,py,Nx,Ny,nof);

% for i=1:Nx
%     for j=1:Ny
%         for k=1:nof
%             Image_patch(:,:,i,j,k) = Image(i:i+px-1,j:j+py-1,k);
%         end
%     end
% end
Image_patch = reshape(Image,[px,Nx,sy,nof]);
Image_patch = permute(Image_patch,[1,3,2,4]);
Image_patch = Image_patch(:,:,:);
Image_patch = reshape(Image_patch,[px,py,Ny,Nx*nof]);
Image_patch = Image_patch(:,:,:);