function [] = piezo_find_calls(directory)
    %Wrapper function for piezo_find_calls_logger
    mergeThresh = 5e-3; % in s
    FS_env = 1000; % Sample frequency of the envelope

    
    files = dir(directory);
    
%     allCallTimes = cell(1, length(files));
%     index = 1;
%     for i = 1:length(files)
%        file = files(i);
%        if contains(file.name, "logger")
%            filepath = strcat(directory, "/", file.name);
%            allCallTimes{index} = piezo_find_calls_logger(filepath);
%            index = index + 1;
%        end
%     end
%     allCallTimes = allCallTimes(1: index - 1);

    %handle any merges    
    a = [1,500,1500,1550];
    b = [2000,2200,2300,2400];
    c = [499,3500];
    test2 = {c, a, b};
    
    allCallTimes = test2;
    answer = zeros(1, sum(cellfun(@length, allCallTimes)));
    i = 1;
    while ~all(cellfun(@isempty, allCallTimes))
        [~, index] = find_first_start(allCallTimes);
        logger = allCallTimes{index};
        allCallTimes{index} = logger(3:end);
        start = logger(1);
        stop = logger(2);
        if i > 1
            prevStop = answer(i - 1);
            if start - prevStop <= mergeThresh * FS_env
               if stop > prevStop 
                    answer(i - 1) = stop;
               end
            else    
                answer(i:i+1) = [start, stop];
                i = i + 2;
            end
        else
            answer(i:i+1) = [start, stop];
            i = i + 2;
        end
        disp(answer)
        disp("-----------------------------------")
    end
    answer = answer(1 : i - 1);
    
end

function [minimum, index] = find_first_start(nested_vectors)
    starts = Inf(1, length(nested_vectors));
    for i = 1:length(nested_vectors)
       tmp = nested_vectors{i};
       if ~isempty(tmp)
            starts(i) = tmp(1);
       end
    end
    [minimum, index] = min(starts);
end

