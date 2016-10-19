function [up down autotime] = calc_crosscorr(filename, combination1, combination2)
% [up down autotime] = calc_crosscorr(filename)
% This is a helper for part_crosscorrelation that is invoked for each data
% package that is generated there

    data = load(filename);

    % We get a matrix, that is composed of column vectors that
    % represent filters of the y data with regards to the arrive
    % time of the photon streams
    %
    % The combinations are indices in this matrix
    for k=1:size(combination1,2)
        filterset1(:, k) = data.filter_matrix(:,combination1(k));
    end
    for k=1:size(combination2,2)
        filterset2(:, k) = data.filter_matrix(:,combination2(k));
    end


    if sum(filterset1(:)) == 0 || sum(filterset2(:)) == 0
        up = zeros(data.Ncasc * data.Nsub, size(filterset1,2), size(filterset2,2));
        down = up;
        autotime = zeros(data.Ncasc * data.Nsub, 1);
        return
    end
    [up, autotime] = tttr2xfcs_cl(data.z, filterset1(:,:), filterset2(:,:), data.Ncasc, data.Nsub);
    % We do a correlation pretending the timescale was reversed - to work
    % around a limitation of the correlation routine (which sorts the y
    % part ascending), we invert by substracting the max(y) value and
    % taking the absolute value -- as the array would the be in the wrong
    % order we invert it and the filter arrays also
    if isfield(data, 'calc_reverse') && data.calc_reverse > 0
        [down] = tttr2xfcs_cl(abs(data.z(end:-1:1)-max(data.z)), filterset1(end:-1:1,:), filterset2(end:-1:1,:), data.Ncasc, data.Nsub);
    else
        [down] = zeros(size(up));
    end
