function [SP_daniel] = Daniel(SP, window)
    SP_daniel = zeros(1, length(SP)); 
    for i=floor(window/2)+1:1:length(SP)-floor(window/2)
        SP_daniel(i) = sum(SP(-floor(window/2)+i:floor(window/2)+i))/window;
    end
end

