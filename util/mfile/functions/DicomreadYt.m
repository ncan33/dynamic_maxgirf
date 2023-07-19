function image = DicomreadYt()

file_name = dir('*.dcm');

for i = 1:size(file_name)
    if ~contains(file_name(i).name(1),'.')
        image(:,:,i) = dicomread(file_name(i).name);
    end
end

idx = sum(sum(image))==0;
image(:,:,idx) = [];

return