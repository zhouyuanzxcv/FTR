function handles = jitter_overlap(handles, b_spacing_multiplier, thresh)
if nargin < 3
    thresh = [0.4, 0.4];
end

size_j = [b_spacing_multiplier, b_spacing_multiplier / 2];

for loopnum = 1:2
    
    for i = 1:length(handles)
        for j = i+1:length(handles)
            h_1 = handles(i);
            h_2 = handles(j);
            h_1_y_interval = get_y_interval(h_1, thresh(2));
            h_2_y_interval = get_y_interval(h_2, thresh(2));
            if intervals_overlap(h_1_y_interval, h_2_y_interval) ...
                    && abs(h_1.XData - h_2.XData) <= thresh(1)
                h_1.XData = h_1.XData - size_j(loopnum);
                h_2.XData = h_2.XData + size_j(loopnum);
            end
        end
    end
end
end

function h_1_y_interval = get_y_interval(h_1, thresh)
h_1_y_interval = [h_1.YData - h_1.YNegativeDelta - thresh/2, ...
                h_1.YData + h_1.YPositiveDelta + thresh/2];
end

function overlap = intervals_overlap(A, B)
% A and B are N×2 matrices where each row is [start, end]
% Returns logical array indicating overlap for each pair

overlap = A(:,1) <= B(:,2) & B(:,1) <= A(:,2);
end

