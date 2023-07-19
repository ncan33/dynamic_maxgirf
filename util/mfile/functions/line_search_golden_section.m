function step = line_search_golden_section(old,update,Data,para)

%step_start = para.Recon.step_size(end)*1.3;
%step_start = 1;
max_try = 10;
golden_ratio = 2/(1+sqrt(5));
%step_old = step_start;
%step_new = step_start;

step_1 = 0;
step_2 = para.Recon.step_size(end);


%cost_1 = para.Cost.totalCost(end);
cost_1 = cost_STCR(old,Data,para);
cost_2 = cost_STCR(old+step_2*update,Data,para);

step_all = zeros(3,2,class(update));
%step_all(:,3) = 1:max_try;
step_all(1,1) = step_1;
step_all(1,2) = cost_1;
step_all(2,1) = step_2;
step_all(2,2) = cost_2;

if cost_2>cost_1
    step_3 = abs(step_2-step_1)*(1-golden_ratio)+min(step_1,step_2);
else
    step_3 = abs(step_2-step_1)*golden_ratio+max(step_2,step_1);
end
step_all(3,1) = step_3;
cost_3 = cost_STCR(old+step_3*update,Data,para);
step_all(3,2) = cost_3;
% now we have 3 valuse, we can start to do golden ratio search!

for i=1:max_try-3
    
    step_sort_cost = sortrows(step_all,2);
    step_sort_step = sortrows(step_all,1);
    %step_sort_cost(step_sort_cost(:,2)==0,:) = [];
    step_min = step_sort_cost(1,1);
    [x,~] = find(step_sort_step==step_min);
    if x==1
        step_1 = step_min;
        step_2 = step_sort_step(2,1);
        step_all(i+3,1) = abs(step_2-step_1)*(1-golden_ratio)+min(step_1,step_2);
        step_all(i+3,2) = cost_STCR(old+step_all(i+3,1)*update,Data,para);
        continue
    end
    if x==i+2
        step_1 = step_min;
        step_2 = step_sort_step(end-1,1);
        step_all(i+3,1) = abs(step_2-step_1)*golden_ratio+max(step_2,step_1);
        step_all(i+3,2) = cost_STCR(old+step_all(i+3,1)*update,Data,para);
        continue
    end
    %step_1 = step_min;
    step_left = step_sort_step(x-1,1);
    step_right = step_sort_step(x+1,1);
    if step_min-step_left > step_right - step_min
        step_all(i+3,1) = abs(step_min-step_left)*(1-golden_ratio)+step_left;
    else
        step_all(i+3,1) = abs(step_right-step_min)*(1-golden_ratio)+step_min;
    end
    step_all(i+3,2) = cost_STCR(old+step_all(i+3,1)*update,Data,para);
    %step_2 = step_sort_cost(2,1);
    %step_3 = abs(step_2-step_1)*(1-golden_ratio)+min(step_1,step_2);
    %step_all(i+3,1) = step_3;

end

step_sort_cost = sortrows(step_all,2);
step = gather(step_sort_cost(1,1));
if step == 0
    keyboard
end
cost = gather(step_sort_cost(1,2));
fprintf(['Step = ' num2str(step) '...'])
fprintf(['Cost = ' num2str(cost) '...'])

function total_cost = cost_STCR(guess,Data,para)

%fidelity= compute_fidelity_yt_new(guess,Data,para);
fidelity_update = bsxfun(@times,guess,Data.sens_map);
fidelity_update = fft2(fidelity_update,para.Recon.kSpace_size(1),para.Recon.kSpace_size(2));
        
fidelity_update_temp(:,:,:,:,1,1,1) = sum(fidelity_update,5);
fidelity_update_temp(:,:,:,:,1,1,2) = fidelity_update(:,:,:,:,1,:,:) + exp(1i*2*pi/3)*fidelity_update(:,:,:,:,2,:,:) + exp(1i*4*pi/3)*fidelity_update(:,:,:,:,3,:,:);
fidelity_update_temp(:,:,:,:,1,1,3) = fidelity_update(:,:,:,:,1,:,:) + exp(1i*4*pi/3)*fidelity_update(:,:,:,:,2,:,:) + exp(1i*2*pi/3)*fidelity_update(:,:,:,:,3,:,:);
fidelity_update_temp = bsxfun(@times,fidelity_update_temp,Data.mask);
        
fidelity_update_temp = Data.kSpace - fidelity_update_temp;
fidelity_cost = sum(abs(fidelity_update_temp(:)).^2);
%fidelity_cost = sum(abs(fidelity(:)).^2);

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