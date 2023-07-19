function step = line_search_bins(old,update,Data,para)

step_start = para.Recon.step_size(end)*1.3;%magic number
%step_start = 2;
%step_start = para.Recon.step_size(1);
tau = 0.8;
max_try = 15;
step = step_start;

cost_old = para.Cost.totalCost(end);

for i=1:max_try
    
    new = old + step*update;
    fidelity_new = compute_fidelity_for_line_search_yt(new,Data,para);
    switch para.Recon.type
        case '3D'
            cost_new = Cost_STCR_3D(fidelity_new,new,para.Recon.weight_sTV,para.Recon.weight_tTV,para.Recon.weight_sliceTV);
        otherwise
            cost_new = Cost_STCR_bins(fidelity_new,new,para);
    end
    
    if cost_new > cost_old
        step = step*tau;
    else
%         fprintf(['Step = ' num2str(step) '...'])
        %fprintf(['Cost = ' num2str(round(cost_new)) '...\n'])
        return
    end

end
% fprintf(['Step = ' num2str(step) '...'])
%fprintf(['Cost = ' num2str(round(cost_new)) '...\n'])
