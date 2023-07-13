% make_patt_bar_opt_starfield_G4
% pattern generator, creates a bar +/- starfield

% INPUTS:
% pattN - pattern number when saving
% objWidth - size of bar
% grtSelect - 0 no starfield, 1 with starfield, -1 against starfield
% starNumber - total number of stars, randomly distributed
% starGS - starfield intensity
% bckGS - background intensity

% 01/11/2022 - MC created

function  make_patt_bar_starfield_pablo_G4(pattN, objWidth, starSelect, starDist, starGS, bckGS)

%% set meta data

userSettings

pattern.x_num = 192 * starDist;
pattern.y_num = 1;

pattern.num_panels = NumofRows*NumofColumns;
pattern.gs_val = 4; %1 or 4
pattern.stretch = zeros(pattern.x_num, pattern.y_num); %match frames
frameN = 16*NumofRows; %px per row
frameM = 16*NumofColumns; %px per col

%initialize pattern data
Pats = zeros(frameN, frameM, pattern.x_num, pattern.y_num);

%set object parameters
objGS = 1;
objPolar = 'dark';


%% generate pattern data

%initialize images
bckImage = ones(frameN, frameM) * bckGS;

barImage = bckImage;
barImage(:,1:objWidth) = objGS;

starImage = bckImage';
tmp = [1:starDist:192] + frameM*[0:starDist:frameN]';
starImage(mod(tmp(:),frameM*frameN)) = starGS;
tmp = circshift(circshift(starImage,ceil(starDist/2),2),ceil(starDist/2),1);
starImage = starImage + tmp;
starImage = starImage';

% store N=starDist frames at each bar position, 1 for each rotation of the
% starField until the pattern repeats
for x = 1:pattern.x_num
    Pats(:,:,x,1) = circshift(starImage,x-1,2);
    Pats(:,mod(ceil(x/starDist)+(1:objWidth)-1,frameM)+1,x,1) = objGS;
end

%save lookup table
patlookup.fullname = [num2str(objWidth) 'px_' num2str(objGS) 'gs' objPolar 'bar' '_with_' num2str(starGS) 'gs_' 'starfield_' num2str(bckGS) 'gsbck'];
patlookup.name = [num2str(objWidth) 'px_' objPolar 'bar' '_with_' 'starfield'];
patlookup.size = num2str(objWidth);
patlookup.object = 'bar + starfield';
patlookup.objgs = objGS;
patlookup.stargs = starGS;
patlookup.bckgs = bckGS;


%store pattern data
pattern.Pats = Pats;
%get the vector data
pattern.data = make_pattern_vector_g4(pattern);


%% save pattern data

%set and save pattern data
pattName = [sprintf('%04d', pattN) '_' patlookup.name];
matFileName = fullfile([exp_path, '\Patterns'], [pattName, '.mat']);
save(matFileName, 'pattern');

%save pattern lookup table
pattLookUp = ['patt_lookup_' sprintf('%04d', pattN)];
matFileName = fullfile(pattern_path, [pattLookUp, '.mat']);
save(matFileName, 'patlookup');

%save pat file
patFileName = fullfile([exp_path, '\Patterns'], ['pat', sprintf('%04d', pattN), '.pat']);
fileID = fopen(patFileName,'w');
fwrite(fileID, pattern.data);
fclose(fileID);
end

