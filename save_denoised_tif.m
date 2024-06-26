all_files = dir('Z:/pablo/stacks/**/*dcdn*.tif');

for i = 1:length(all_files)
    if isempty(dir([all_files(i).folder,'\registration_001\denoised_regProduct.mat']))
    fprintf('%s\n',all_files(i).folder)

    tmp = dir([all_files(i).folder,'\*metadatanew*.mat']);
    load([tmp.folder,'\',tmp.name])

    pth_raw_tif = [all_files(i).folder,'\',all_files(i).name];
    
    sz = [md.ypix md.xpix md.numslice md.numvol];
    size_z_read_from_raw = md.numslice_withflyback; %raw tif includes flyback
    size_z_read_from = sz(3);
    size_t_read_from = sz(4);
    inds_z_read_from = 1:size_z_read_from; %can choose to not read the flyback frames here
    inds_t_read_from = 1:size_t_read_from;
    size_read_to = [length(inds_t_read_from), length(inds_z_read_from) sz(1) sz(2)]; %read the way it was written for speed, permute within cx_read_tif_tzyx
    
    regProduct = nan(sz);
    
    tr = Tiff(pth_raw_tif);
    
    for t = 1:size_t_read_from
        for z = 1:size_z_read_from
            tr.setDirectory(double(z + size_z_read_from*(t-1)))
            regProduct(:,:,z,t) = tr.read();
        end
    end
    
    save([all_files(i).folder,'\registration_001\denoised_regProduct.mat'],'regProduct','-v7.3')
    end
end