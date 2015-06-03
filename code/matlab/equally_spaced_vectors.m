function [idx]=equally_spaced_vectors(uvnodell,uvnode,nodexy,trinodes,region,spacing)





[tld tli]=sort((uvnodell(:,2)-region(4)).^2 + (uvnodell(:,1)-region(1)).^2);
[trd tri]=sort((uvnodell(:,2)-region(4)).^2 + (uvnodell(:,1)-region(2)).^2);
[bld bli]=sort((uvnodell(:,2)-region(3)).^2 + (uvnodell(:,1)-region(1)).^2);


xmultiplier=floor(abs(uvnode(tli(1),1)-uvnode(tri(1),1))/spacing)+11
ymultiplier=floor(abs(uvnode(tli(1),2)-uvnode(bli(1),2))/spacing)+14

	
XI=(uvnode(tli(1),1):spacing:(uvnode(tli(1),1)+(xmultiplier*spacing)));
YI=(uvnode(tli(1),2):-spacing:(uvnode(tli(1),2)-(ymultiplier*spacing)));
[XI,YI]=meshgrid(XI,YI); 
	

[a b]=size(XI);
idx=zeros(a*b,1);
	


cnt=1;
for i=1:length(trinodes)
	if max(max(inpolygon(XI,YI,nodexy(trinodes(i,:),1),nodexy(trinodes(i,:),2))))==1
		idx(cnt)=i;
		cnt=cnt+1;
	end
end
	
idx(idx==0) = [];


