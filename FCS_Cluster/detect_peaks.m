function bounds = detect_peaks(x_unsorted, y_unsorted, count)
% function bounds = detect_peaks(x_unsorted, y_unsorted, count)
%
% This algorithm takes two column vectors and one scalar. The
% vectors describe a histogramm of the values (x) with respect to
% their count (y). In this histogram count peaks are searchend and
% reported as first column of the result (bounds). The second column
% reports the minimum behind the corresponding peak, but prior
% to the next peak.
%
% Potential problems:
% - to small peak distance (has to be at least max(x(y>0)) / (1.5 * count))

    % BEWARE -- until the very end -- bounds are noted as indices into x
    % at the end bounds are converted to the "real" x data.

    %x = x_unsorted;
    y = y_unsorted;
    % Make sure we have the values sorted ascending in regards to X
    [x, x_ind] = sort(x_unsorted);
    %y = zeros(size(y_unsorted));
    for i=1:(size(x_ind, 1) * size(x_ind, 2))
        y(i) = y_unsorted(x_ind(i));
    end
    
    % Initalize bounds matrix
    bounds = zeros(count,2);
    bounds(:,1) = -Inf;
    bounds(:,2) = Inf;
    
    found = 0;
    
    % Define the minimum distance between peaks
    min_distance = max(x(y>0)) / (1.5 * count);
    
    % make a copy of y, that we can play with
    search_space = y;
    
    while sum(search_space) > 0 && found < count
        % Find the maxima of the rest of the search space
        candidates = find(search_space == max(search_space));
        % Delete these from search space => next round we find the next
        % lower ones
        search_space(candidates) = zeros;
        % Sort descending - get the most right peak first
        candidates = sort(candidates, 1, 'descend');
        % Check, that the candiates lie outside the already found peaks
        % (min_distance is the space that is kept free around a found peak
        for i=1:size(candidates)
            candidate = candidates(i);
            ok = 1;
            for j=1:count
                if bounds(j,1) ~= -Inf && ...
                   x(candidate) < (x(bounds(j,1)) + min_distance)  && ...
                   x(candidate) > (x(bounds(j,1)) - min_distance)
                   ok = 0;
                   break;
                end
            end
            % If this candidate is ok - see above - add the index to the
            % bounds we found
            if ok
                found = found + 1;
                bounds(found, 1) = candidate;
                break;
            end
        end
    end
    
    bounds(:,1) = sort(bounds(:,1));
    
    for i=1:count
        lower_bound = bounds(i, 1);
        if i==count
            upper_bound = find(x == max(x));
        else
            upper_bound = bounds(i+1, 1);
        end
        upper_bound = upper_bound(1);
        min_indices = find(y == min(y(lower_bound:upper_bound)));
        for j=upper_bound:-1:lower_bound
            if sum(min_indices == j) > 0
                bounds(i,2) = j;
            end
        end
    end
    
    % Convert from indices into x vector into real data
    for i=1:(size(bounds,1) * size(bounds,2))
        bounds(i) = x(bounds(i));
    end