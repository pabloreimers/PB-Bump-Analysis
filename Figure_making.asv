%% Plot the r2 of linear fits for forward, rotation, and joint velocity fits
tmp_idx = lpsp_idx & include_idx;
tmp = T(tmp_idx,:);


figure(1); clf
swarmchart([1,2,3].*[ones(height(tmp),3)],tmp{:,end-2:end},'filled','xjitterwidth',.25,'MarkerFaceColor',.5*[1,1,1])
xticks([1,2,3])
xticklabels({'Forward','Rotational','Joint'})
ylabel('R^2')
title('Linear Fits dF/F and Velocity')

