function step = line_search_reorder(old,update,Data,para)

step_start = para.Recon.step_size(end)*1.3;%magic number
%step_start = 2;
%step_start = para.Recon.step_size(1);
tau = 0.8;
max_try = 15;
step = step_start;

cost_old = para.Cost.totalCost(end);
flag = 0;

for i=1:max_try
%     fprintf(['Iter = ' num2str(i) '...\n'])
    
    new = old + step*update;
    fidelity_new = compute_fidelity_for_line_search_yt(new,Data,para);
    switch para.Recon.type
        case {'3D','3D SOS','3D SOS server'}
            cost_new = Cost_STCR_3D(fidelity_new,new,para.Recon.weight_sTV,para.Recon.weight_tTV,para.Recon.weight_sliceTV);
        otherwise
            cost_new = Cost_STCR(fidelity_new,new(Data.order),para.Recon.weight_sTV,para.Recon.weight_tTV);
    end
    
    if cost_new > cost_old && flag == 0
        step = step*tau;
    elseif cost_new < cost_old 
        step = step*1.3;
        cost_old = cost_new;
        flag = 1;
    elseif cost_new > cost_old && flag == 1
        step = step/1.3;
%         fprintf(['Step = ' num2str(step) '...\n'])
%         fprintf(['Cost = ' num2str(round(cost_old)) '...\n'])
        return
    end
end
% fprintf(['Step = ' num2str(step) '...\n'])
% fprintf(['Cost = ' num2str(round(cost_new)) '...\n'])



function total_cost = cost_TCR(guess,Data,para)

fidelity= compute_fidelity_yt_new(guess,Data,para);
fidelity_cost = sum(abs(fidelity(:)).^2);        

if para.weight_tTV~=0
    temporal_cost = diff(guess,1,3);
    temporal_cost = para.weight_tTV*sum(abs(temporal_cost(:)));
else
    temporal_cost = 0;
end
if para.weight_sTV~=0
    spatial_cost_x = diff(guess,1,1);
    spatial_cost_y = diff(guess,1,2);
    spatial_cost_x(end+1,:,:,:,:,:) = 0;
    spatial_cost_y(:,end+1,:,:,:,:) = 0;
    spatial_cost = sqrt(abs(spatial_cost_x).^2+abs(spatial_cost_y).^2);
    spatial_cost = para.weight_sTV*sum(spatial_cost(:));
else
    spatial_cost = 0;
end
total_cost = fidelity_cost+temporal_cost+spatial_cost;
total_cost = total_cost/numel(guess);