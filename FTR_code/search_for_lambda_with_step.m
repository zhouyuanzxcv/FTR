function [lambda_best, extra_best, obj_vals, best_ind] = search_for_lambda_with_step(...
    params, lambdas, start_lambda, eval_fun)


obj_vals = NaN(1, length(lambdas));
extras = cell(1, length(lambdas));


pos_curr = find(lambdas == start_lambda); % start at 2l which should be less than l1_max

[pos_curr, obj_vals, extras] = search_for_lambda_loop(params, ...
    lambdas, pos_curr, obj_vals, extras, eval_fun);

best_ind = pos_curr;

lambda_best = lambdas(pos_curr);

extra_best = extras{pos_curr};



end

function [pos_curr, obj_vals, extras] = search_for_lambda_loop(params, ...
    lambdas, pos_curr, obj_vals, extras, eval_fun)

while 1
    [obj_val_curr, obj_vals, extras] = get_obj_val_by_lambda(params, ...
        lambdas, pos_curr, obj_vals, extras, eval_fun);

    pos_right = min(pos_curr + 1, length(lambdas));
    [obj_val_right, obj_vals, extras] = get_obj_val_by_lambda(params, ...
        lambdas, pos_right, obj_vals, extras, eval_fun);

    pos_left = max(pos_curr - 1, 1);
    [obj_val_left, obj_vals, extras] = get_obj_val_by_lambda(params, ...
        lambdas, pos_left, obj_vals, extras, eval_fun);

    if obj_val_curr <= obj_val_left && obj_val_curr <= obj_val_right
        break;
    else
        if obj_val_left < obj_val_right
            pos_curr = pos_left;
        else
            pos_curr = pos_right;
        end
    end

end

end

function [obj_val, obj_vals, extras] = get_obj_val_by_lambda(params, ...
    lambdas, pos_curr, obj_vals, extras, eval_fun)
if isnan(obj_vals(pos_curr)) % current value not evaluated yet
    [obj_val, extra] = eval_fun(params, lambdas(pos_curr));
    obj_vals(pos_curr) = obj_val;
    extras{pos_curr} = extra;
else
    obj_val = obj_vals(pos_curr);
end

end

