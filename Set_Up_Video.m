function v = Set_Up_Video(videoname)    

    % Prepare video writer
    todayStr = datestr(now, 'yyyy-mm-dd');
    saveFolder = fullfile('C:\Users\jewellan\Documents\MATLAB\ORACLE\Results', todayStr);
    
    if ~exist(saveFolder, 'dir')
        mkdir(saveFolder);
    end
    
    filePath = fullfile(saveFolder, videoname);
    v = VideoWriter(filePath, 'MPEG-4');
    v.FrameRate = 10; % adjust playback speed

end