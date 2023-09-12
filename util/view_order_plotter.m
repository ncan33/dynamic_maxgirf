function view_order_plotter(kspace_info)
    view_order = kspace_info.viewOrder;
    new_view_order = kspace_info.viewOrder;
    counter = 0;
    for i = 1:size(view_order, 2)
        if view_order(i) > 10
            new_view_order(i) = view_order(i) - 10;
            counter = counter + 1;
        end
    end
    
    %plot(new_view_order)
    plot(new_view_order(1:30))
end