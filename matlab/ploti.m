function dcm_obj = ploti(a)
% Plots graph and sets up a custom data tip update function
fig = figure('DeleteFcn','doc datacursormode');
scatter(a(:,1), a(:,2));
dcm_obj = datacursormode(fig);
set(dcm_obj,'UpdateFcn',{@myupdatefcn});

function txt = myupdatefcn(~,event_obj)
pos = get(event_obj,'Position');
I = get(event_obj, 'DataIndex');
txt = {['X: ',num2str(pos(1))],...
      ['Y: ',num2str(pos(2))],...
      ['I: ',num2str(I)]};