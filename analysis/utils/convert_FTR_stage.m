function stage2 = convert_FTR_stage(subtype_stage, traj, mode)
reparam_mode = 'L2';
nsubtype = size(traj, 3);
num_int = size(traj,2);
ds = zeros(num_int, nsubtype);

for k = 1:nsubtype
    f0 = traj(:,:,k)'; % f0 already starts at 0

    if strcmp(reparam_mode, 'L2')
        df = f0(2:end,:) - f0(1:end-1,:);
        d = sqrt(sum(df.^2, 2));
    elseif strcmp(reparam_mode, 'L1')
        df = f0(2:end,:) - f0(1:end-1,:);
        d = sum(abs(df), 2);
    end

    d = [0;d];   % point 1 to the origin is 0
    d = cumsum(d);
    ds(:,k) = d;
end

subtype = subtype_stage(:,1);
stage = subtype_stage(:,2);
stage2 = stage;

normalized_stage = linspace(0, 1, num_int);

for k = 1:nsubtype
    stage_k = stage(subtype == k);
    if strcmp(mode, 'normalized2absolute')
        stage_idx = stage_k * (num_int - 1) + 1;
        stage_k_new = ds(round(stage_idx), k);
    elseif strcmp(mode, 'absolute2normalized')
        diff = repmat(stage_k, 1, num_int) - repmat(ds(:,k)', length(stage_k), 1);
        [~,min_idx] = min(abs(diff), [], 2);
        stage_k_new = normalized_stage(min_idx);
    end
    stage2(subtype == k) = stage_k_new;
end


end
