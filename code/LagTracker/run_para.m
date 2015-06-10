%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 
% Return: 
%
%   Run code in parallel and save particle postions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_times=12

matlabpool open 4
	parfor starti=1:start_times
		starti
        pause(starti)
		savelag{starti}=main(starti);
	end
savepath = ['/array/home/119865c/karsten/lagtracker/savedir/gp_aug01_2013.mat']
save(savepath,'savelag');

matlabpool close
