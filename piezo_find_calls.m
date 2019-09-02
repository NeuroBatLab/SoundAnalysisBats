function [] = piezo_find_calls(directory)
    % Wrapper function for piezo_find_calls_logger. Takes in a directory
    % (should be named audiologgers) and calls piezo_find_calls_logger.m
    % on each logger in the directory. Then takes the output of those calls
    % and merges the call times for any calls that are less than
    % mergeThresh. Note this latter behavior is subject to change.
    
    mergeThresh = 5e-3; % minimum length between two calls in s
    FS_env = 1000; % Sample frequency of the envelope

    files = dir(directory); % get all of the files in the directory

    % Go through all of the files in the directory and call
    % piezo_find_calls_logger.m on them if they are a logger data file. 
    allCallTimes = cell(1, length(files));
    index = 1;
    for ii = 1:length(files)
       file = files(ii);
       if contains(file.name, "logger")
           filepath = strcat(directory, '/', file.name, '/', 'extracted_data', '/');
           allCallTimes{index} = piezo_find_calls_logger(filepath);
           index = index + 1;
       end
    end
    allCallTimes = allCallTimes(1: index - 1);

    % test data set
%     a = [1, 500];
%     b = [1500, 1550];
%     logger1 = {a, b};
%     c = [2000, 2200];
%     d = [2300, 2400];
%     logger2 = {c, d};
%     e = [499, 3500];
%     logger3 = {e};
%     test2 = {logger3, logger2, logger1};
%     
%     allCallTimes = test2;

    % Handle the case of any calls that should be merged together. Note
    % this behavior is subject to change as it seems a bit weird to do it
    % like this. 
    answer = zeros(1, sum(cellfun(@length, allCallTimes)));
    ii = 1;
    while ~all(cellfun(@isempty, allCallTimes))
        [~, index] = find_first_start(allCallTimes);
        logger = allCallTimes{index};
        if iscell(logger)
        logger = logger{1};
        end
        allCallTimes{index} = [];
        start = logger(1);
        stop = logger(2);
        if ii > 1
            prevStop = answer(ii - 1);
            if start - prevStop <= mergeThresh * FS_env
               if stop > prevStop 
                    answer(ii - 1) = stop;
               end
            else    
                answer(ii:ii+1) = [start, stop];
                ii = ii + 2;
            end
        else
            answer(ii:ii+1) = [start, stop];
            ii = ii + 2;
        end
        disp(answer)
    end
    answer = answer(1 : ii - 1);
    
end

function [minimum, index] = find_first_start(nested_vectors)
    % Takes in a cell of callTimes (itself a cell of start/stop times) 
    % NESTED_VECTORS and returns the index and minimum value of the logger 
    % with the first start time.

    starts = Inf(1, length(nested_vectors));
    for ii = 1:length(nested_vectors)
       logger = nested_vectors{ii};
       if iscell(logger)
        logger = logger{1};
       end
       if ~isempty(logger)
            starts(ii) = logger(1);
       end
    end
    [minimum, index] = min(starts);
end

