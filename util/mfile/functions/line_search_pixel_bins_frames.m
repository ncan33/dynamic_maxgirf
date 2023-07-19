function step = line_search_pixel_bins_frames(old,update,Data,para)

nof = size(old,3);
step_start = ones(1,1,nof);
step_flag = logical(step_start);

tau = 0.8;
max_try = 15;
step = step_start;

fidelity_new = compute_fidelity_yt_new(old,Data,para);
cost_old = Cost_STCR_pixel_bins_frames(fidelity_new,old,Data,para);

for i=1:max_try
    
    new = old + step.*update;
    fidelity_new = compute_fidelity_yt_new(new,Data,para);
    cost_new = Cost_STCR_pixel_bins_frames(fidelity_new,new,Data,para);
    
    for nf = 1:nof
        if step_flag(1,1,nf) && cost_new(1,1,nf) > cost_old(1,1,nf)
            step(1,1,nf) = step(1,1,nf)*tau;
        else
            step_flag(1,1,nf) = false;
        end
    end
    
    if sum(step_flag) == 0
        figure,plot(squeeze(step))
        return
    end

end
step(step_flag) = 0;
%fprintf(['Step = ' num2str(step) '...'])
%fprintf(['Cost = ' num2str(round(cost_new)) '...\n'])



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