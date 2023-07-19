function Make_Video_All(folder,save_dir)
mkdir(save_dir)
if folder(end) ~= '/'
    folder(end+1) = '/';
end
if save_dir(end) ~= '/'
    save_dir(end+1) = '/';
end
sp = find(folder == '/');
Patient = folder(sp(end-1)+1:sp(end)-1);
mat_files = dir([folder,'ReconData/mat_files/*.mat']);

N = length(mat_files);

save_dir = [save_dir,Patient];

Orin_lookup = 'MID_Orin Look Up';

for i=1:N
    try
        load([mat_files(i).folder,'/',mat_files(i).name],'Image')

        MID = strfind(mat_files(i).name,'MID');
        MID = mat_files(i).name(MID+3:MID+5);
        if strfind(Orin_lookup,MID)
            Orin = str2num(Orin_lookup(strfind(Orin_lookup,MID)+4));
        else
            Orin_lookup(end+2:end+4) = MID;
            Orin = 0;
        end

        if length(size(Image)) == 4
            Orin = Make_3D_video(Image,[save_dir,'_',mat_files(i).name(1:end-18),'.avi'],Orin);
            Orin_lookup(strfind(Orin_lookup,MID)+4) = num2str(Orin);
        elseif length(size(Image)) == 3
            Orin = Make_2D_video(Image,[save_dir,'_',mat_files(i).name(1:end-18),'.avi'],Orin);
            Orin_lookup(strfind(Orin_lookup,MID)+4) = num2str(Orin);
        else
            error('Image is not 2D or 3D.')
        end
    catch
    end
    close all
end
end

function Orin = Make_3D_video(Image,file_name,Orin)
    
    [sx,sy,sz,nf] = size(Image);
    
    sx_cut = round(sx/4);
    Image(1:sx_cut,:,:,:) = [];
    Image(end-sx_cut+1:end,:,:,:) = [];
    sy_cut = round(sy/4);
    Image(:,1:sy_cut,:,:) = [];
    Image(:,end-sy_cut+1:end,:,:) = [];
    
    Image = abs(Image);
    
    if mod(sz,2) ~= 0
        Image(:,:,end+1,:) = zeros(sx,sy,1,nf);
    end
    
    if Orin == 0
        Orin = orintation_detection(sum(Image(:,:,round(sz/2),:),4));
    end
    
    Image = orintate_image(Image,Orin);
    
    [sx,sy,sz,nf] = size(Image);
    Image = reshape(Image,[sx,sy,2,sz/2,nf]);
    Image = permute(Image,[1 3 2 4 5]);
    Image = reshape(Image,[sx*2,sy*sz/2,nf]);
    
    v = VideoWriter(file_name);
    v.FrameRate = 10;
    v.Quality = 100;
    open(v)
    
    Video = figure;
    for n = 1:nf
        imagesc(Image(:,:,n))
        colormap gray
        axis image
        axis off
        brighten(0.4)
        drawnow
        image = getframe;
        writeVideo(v,image)
    end
    close(v)
    close(Video)
end

function Orin = Make_2D_video(Image,file_name,Orin)
        
    [sx,sy,nf] = size(Image);
    
    sx_cut = round(sx/4);
    Image(1:sx_cut,:,:,:) = [];
    Image(end-sx_cut+1:end,:,:,:) = [];
    sy_cut = round(sy/4);
    Image(:,1:sy_cut,:,:) = [];
    Image(:,end-sy_cut+1:end,:,:) = [];
    
    Image = abs(Image);
    
    if Orin == 0
        Orin = orintation_detection(sum(Image,3));
    end
    
    Image = orintate_image(Image,Orin);
    
    v = VideoWriter(file_name);
    v.FrameRate = 10;
    v.Quality = 100;
    open(v)
    
    Video = figure;
    for n = 1:nf
        imagesc(Image(:,:,n))
        colormap gray
        axis image
        axis off
        brighten(0.4)
        drawnow
        image = getframe;
        writeVideo(v,image)
    end
    close(v)
    close(Video)
end