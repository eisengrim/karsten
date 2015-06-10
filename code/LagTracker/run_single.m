%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 
% Return: 
%
%   Run code on single core and save particle postions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% savepath = ['/array/home/119865c/karsten/lagtracker/savedir/gp_aug01_2013_s' num2str(starti) '.mat']
start_times=12;

%single cpu 
	for starti=0:3:21
		starti
        clear savelag
        savelag=main(starti);
        savelag=struct(savelag)
        
        savepath = ['/array/home/119865c/karsten/lagtracker/savedir/gp_aug01_2013_s' num2str(starti) '.mat']
        % save(['savedir/kit4_kelp_nodrag/kit4_kelp_tight2_kelpfield_3elements_200x200_1000pp_s' num2str(starti) '.mat'],'savelag','-v7.3');
        save(savepath, 'savelag', '-v7.3');
	end

