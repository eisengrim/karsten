function [region, passage,runstart]=regions(name)



runstart=datenum([0 0 0 0 0 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isequal(name,'pp')
		runstart=datenum([2012 09 12 0 0 0]); %pp start time
		region=[-66.225 -66.195 44.37 44.41];
		passage='Petit Passage';
	end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%gp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isequal(name,'gp')
		runstart=datenum([2012 07 21 0 0 0]); %gp start time
		region=[-66.38 -66.29 44.213 44.32];
		passage='Grand Passage';
	end
	if isequal(name,'gp_tight')
		runstart=datenum([2012 07 21 0 0 0]); %gp start time	
		region=[-66.355 -66.32 44.245 44.2925];
		passage='Grand Passage Tight';
	end
	if isequal(name,'gp3')
		region=[-66.345 -66.33 44.26 44.275];
		passage='Grand Passage Site 3';
	end
	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isequal(name,'dg')
		runstart=datenum([2012 01 25 0 0 0]); %dg start time
		region=[-65.79 -65.73 44.65 44.7];
		passage='Digby Gut';
	end
	if isequal(name,'dg_upper')
		runstart=datenum([2012 01 25 0 0 0]); %dg start time
		region=[-65.794 -65.743 44.66 44.71];
		passage='Digby Gut Upper';
	end	




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mp 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isequal(name,'big_force')	
		region=[-64.455 -64.37 45.353 45.383];%forcesite to cape sharp
		passage='Big Force Site';
	end
	if isequal(name,'force_plus')		
		region=[-64.4435 -64.40 45.355 45.3715];%forcesite plus more to east and south
		passage='Force Site';
	end	
	if isequal(name,'force')		
		region=[-64.4435 -64.42 45.3615 45.3715];%forcesite
		passage='Force Site';
	end	
	if isequal(name,'mp')
		region=[-64.52 -64.3 45.3 45.4];
		passage='Minas Passage';
	end
	if isequal(name,'channel')
		region=[-65 -64 45.1 45.5];
		passage='Minas Channel';
	end
	if isequal(name,'br')
		region=[-64.425 -64.40 45.36 45.375];
		passage='Black Rock';
	end
	if isequal(name,'capesplit')
		region=[-64.515 -64.48 45.33 45.345];
		passage='Cape Split';
	end
	if isequal(name,'sharpbr')
		region=[-64.440 -64.360 45.347 45.375];
		passage='Black Rock and Cape Sharp';
	end	
	if isequal(name,'westbay')
		region=[-64.425  -64.325   45.3547   45.385];
		passage='West Bay';
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%general
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isequal(name,'dn')
		runstart=datenum([2012 07 21 0 0 0]); %gp start time
		region=[-66.4 -65.7 44.22 44.75];
		passage='Digby Neck';
	end
	if isequal(name,'dn_all')
		region=[-66.45 -65.5 44.20 44.75];
		passage='Digby Neck';
	end
	if isequal(name,'maritimes')
		region=[-70 -59 41 49];
		passage='Maritimes';
	end
	if isequal(name,'ns')
		region=[-66.5 -59.6 43.25 47.25];
		passage='Nova Scotia';
	end
	if isequal(name,'abasin')		
		region=[-65.81 -65.43 44.575 44.79];
		passage='Annapolis Basin';
	end
	if isequal(name,'sfmwhole')
		region=[-72   -56    37    47];
		passage='Nova Scotia';
	end
































