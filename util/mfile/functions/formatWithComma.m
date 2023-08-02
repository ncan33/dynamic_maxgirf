function formattedStr = formatWithComma(num)
    % Convert the number to a string representation
    str = num2str(num);
    
    % Determine the position to insert the comma
    commaPos = mod(numel(str) - 1, 3) + 1;
    
    % Insert the comma as a thousand separator
    formattedStr = str;
    while commaPos < numel(formattedStr)
        formattedStr = [formattedStr(1:commaPos), ',', formattedStr(commaPos+1:end)];
        commaPos = commaPos + 4; % Move to the next comma position
    end
end