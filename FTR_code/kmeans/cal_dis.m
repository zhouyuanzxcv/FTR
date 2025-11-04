function [stag,dist,diff] = cal_dis(dat,traj)
    [stag,dist,diff] = cal_dis_1(dat,traj);
end


function [stag,dist,diff] = cal_dis_1(dat,traj)
    num_int = size(traj,2);
    [~,I] = pdist2(traj.',dat,'euclidean','Smallest',1);
    stag = (I-1)/(num_int-1);
    dist = (dat-traj(:,I).').^2;
    diff = dat-traj(:,I).';
end

function [stag,dist,diff] = cal_dis_2(dat,traj)
    npoin = size(traj,2);
    nsamp = size(dat,1);
    max_ep = ceil(log2(npoin));
    head = ones(nsamp,1);
    tail = ones(nsamp,1)*npoin;
    mid = ceil(1/2*(head+tail));
    for ep = 1:max_ep
        comp = sum((dat-traj(:,head).').^2,2) < sum((dat-traj(:,tail).').^2,2);
        tail = mid.*comp + tail.*(1 - comp);
        head = head.*comp + mid.*(1 - comp);
        mid = ceil(1/2*(head+tail));
    end
    dist = (dat-traj(:,mid).').^2;
    diff = dat-traj(:,mid).';
    stag = mid/npoin;
end