function step = line_search_sr(old,update,Data,Data_off,para)

step_start = para.Recon.step_size(end)*1.3;%magic number
%step_start = 2;
%step_start = para.Recon.step_size(1);
tau = 0.8;
max_try = 15;
step = step_start;

cost_old = para.Cost.totalCost(end);
order_back = para.Recon.order_back;
nslice_low_res = para.Recon.nslice_low_res;
siz = para.Recon.siz;

for i=1:max_try
    
    new = old + step*update;
    
    Image_1 = reshape(new(:,:,:,:,1:end-1),[siz(1:4),2,nslice_low_res]);
    Image_1 = permute(sum(Image_1,5),[1,2,3,4,6,5]);
    Image_2 = reshape(new(:,:,:,:,2:end),[siz(1:4),2,nslice_low_res]);
    Image_2 = permute(sum(Image_2,5),[1,2,3,4,6,5]);
    
    for j=1:para.Recon.nset
        fidelity_1(j) = compute_fidelity_for_line_search_yt(Image_1(:,:,:,:,order_back(:,j)),Data{j},para);
        fidelity_2(j) = compute_fidelity_for_line_search_yt(Image_2(:,:,:,:,order_back(:,j)),Data_off{j},para);
    end
    fidelity_new = sos([fidelity_1,fidelity_2])/10^2;
    
%    fidelity_new = compute_fidelity_for_line_search_yt(new,Data,para);
%     switch para.Recon.type
%         case '3D'
%             cost_new = Cost_STCR_3D(fidelity_new,new,para.Recon.weight_sTV,para.Recon.weight_tTV,para.Recon.weight_sliceTV);
%         otherwise
            cost_new = Cost_STCR(fidelity_new,new,para.Recon.weight_sTV,para.Recon.weight_tTV);
%     end
    
    if cost_new > cost_old
        step = step*tau;
    else
        %fprintf(['Step = ' num2str(step) '...\n'])
        %fprintf(['Cost = ' num2str(round(cost_new)) '...\n'])
        return
    end

end
%fprintf(['Step = ' num2str(step) '...\n'])
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