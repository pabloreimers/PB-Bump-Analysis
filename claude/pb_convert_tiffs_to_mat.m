%PB_CONVERT_TIFFS_TO_MAT Batch-convert denoised TIFF movies to .mat files.
%
%   For every 'denoised_srdtrans.tif' found anywhere under ROOT, read it with
%   pb_read_tiff_movie and save the [Y,X,T] array as 'denoised_srdtrans.mat'
%   (variable name: movie) in the same folder as the tiff.
%
%   Set OVERWRITE = true to re-convert trials that already have a .mat.
%
%   Run pb_read_tiff_movie.m must be on the path (it lives in this folder).

root      = 'Z:\pablo\epg_gain_change\';
overwrite = false;   % skip trials that already have a .mat unless this is true

hits = dir(fullfile(root, '**', 'denoised_srdtrans.tif'));
nHits = numel(hits);
fprintf('pb_convert_tiffs_to_mat: found %d tiff(s) under %s\n\n', nHits, root);

nDone = 0; nSkip = 0; nFail = 0;
failed = {};

for k = 1:nHits
    tifPath = fullfile(hits(k).folder, hits(k).name);
    matPath = fullfile(hits(k).folder, 'denoised_srdtrans.mat');

    fprintf('[%d/%d] %s\n', k, nHits, hits(k).folder);

    if ~overwrite && isfile(matPath)
        fprintf('  skip: .mat already exists\n\n');
        nSkip = nSkip + 1;
        continue;
    end

    try
        movie = pb_read_tiff_movie(tifPath); %#ok<NASGU>
        % -v7.3 handles arrays larger than 2 GB (movies can exceed this).
        save(matPath, 'movie', '-v7.3');
        clear movie;
        fprintf('  saved: %s\n\n', matPath);
        nDone = nDone + 1;
    catch ME
        fprintf(2, '  FAILED: %s\n\n', ME.message);
        failed{end+1} = tifPath; %#ok<SAGROW>
        nFail = nFail + 1;
    end
end

fprintf('=== done: %d converted, %d skipped, %d failed (of %d) ===\n', ...
    nDone, nSkip, nFail, nHits);
if nFail > 0
    fprintf(2, 'Failed trials:\n');
    fprintf(2, '  %s\n', failed{:});
end
