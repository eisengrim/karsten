PATH2SIM='/EcoII/acadia_uni/workspace/simulated/FVCOM/dngridCSR/drifter_runs/GP/2013_Aug_01_3D/output/subdomain_GP1_0001.nc'
PATH2OBS='/EcoII/acadia_uni/workspace/observed/GP/Drifter/GP_F_20130801_78_2_001_sE15.mat'

if 1==1
	
	name='gp_d_aug01_13';
	grid='dngridCSR';	
	type='2d';
	resolution='-r600';	
	

    load(PATH2SIM);
    data=load(PATH2OBS);

    region=[min(nodell(:,1)) max(nodell(:,1)) min(nodell(:,2)) max(nodell(:,2))];

    ua=squeeze(u(:,1,:))';
    va=squeeze(v(:,1,:))';
end

speed=sqrt(ua.^2+va.^2);

[idx]=equally_spaced_vectors(uvnodell,uvnode,nodexy,trinodes,region,100);

loc=['~/karsten/data/gp/' grid '_' name '/'];
mkdir(loc);


%plot time-steps
if 1==1
	figure
	for j=71:length(time)
	j
	

	
		%plot speed
		clf
		%axis(region)
	
		%plot_google_map('maptype','satellite')
		patch('Vertices',nodell,'Faces',trinodes,'Cdata',speed(j,:))
		shading flat
		axis(region)
		caxis([0 3])
		colorbar		
		
		hold on
	 		scale=.0005;
			uvxy=uvnodell(idx,:);
    		offset=[ ua(j,idx); va(j,idx) ]';
			arrow(uvxy,uvxy+(offset*scale),'Length',2,'BaseAngle',25,'TipAngle',35,'Width',.6)

 		scale=.0005;
            [b a]=find(data.velocity.vel_time>(time(j)-datenum([0 0 0 0 5 0])) & data.velocity.vel_time<(time(j)+datenum([0 0 0 0 5 0])));
      

            if (isempty(a))
            'empty'
            else
                a=a(1:10:end)
                uvxy_d=[data.velocity.vel_lon(a); data.velocity.vel_lat(a)];
                offset_d=[ data.velocity.u(a); data.velocity.v(a) ];
			    arrow(uvxy_d,uvxy_d+(offset_d*scale),'Length',2,'BaseAngle',25,'TipAngle',35,'Width',.6,'FaceColor',[1 1 1], 'EdgeColor',[1 1 1])
            end
		hold off

	
	
		title(['' datestr( time(j) , 'YYYY mmmm DD - HH:MM:SS' ) ''])
		set(gca,'fontsize',18)	
		plotbox = get_aspectratio([region(1) region(2)],[region(3) region(4)],1);
		set(gca, 'PlotBoxAspectRatio', plotbox)

		name2=['' grid '_' name '_speed_step_' num2str(j,'%03d') ''];
		%print('-dpng',resolution,['' loc '' name2 '.png']);	
		print('-dpng',resolution,['' loc '' name2 '_with_drifter.png']);





	end

end








