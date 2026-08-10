run('..\data scripts\epg_dlight_claude_script.m');
for k = 3:5
    saveas(figure(k), sprintf('fig%d.png',k));
end
