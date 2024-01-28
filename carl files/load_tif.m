[file,path] = uigetfile('.mat','Select metadata');
load([path,'\',file])
[file,path] = uigetfile(path,'Select TIF');
pth_raw_tif = [path,'\',file];

out_datatype = "int16"; %NOTE: SCANIMAGE SAVES INT16 NOT UINT16, CONVERT BELOW
%%
sz = [md.ypix md.xpix md.numslice md.numvol];
size_z_read_from_raw = md.numslice_withflyback; %raw tif includes flyback
size_z_read_from = sz(3);
size_t_read_from = sz(4);
inds_z_read_from = 1:size_z_read_from; %can choose to not read the flyback frames here
inds_t_read_from = 1:size_t_read_from;
size_read_to = [length(inds_t_read_from), length(inds_z_read_from) sz(1) sz(2)]; %read the way it was written for speed, permute within cx_read_tif_tzyx

imgData = nan(sz);

tr = Tiff(pth_raw_tif);

for t = 1:size_t_read_from
    for z = 1:size_z_read_from
        tr.setDirectory(double(z + size_z_read_from*(t-1)))
        imgData(:,:,z,t) = tr.read();
    end
end