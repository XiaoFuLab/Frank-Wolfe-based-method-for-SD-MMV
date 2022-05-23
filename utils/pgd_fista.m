function [x_k, tracking] = pgd_fista(g_fn, p_fn, step_size, init_point, ops) 

    if isfield(ops, 'max_iters'), max_iters = ops.max_iters; else max_iters = 500; end;
    if isfield(ops, 'debug'), DEBUG = ops.debug; else DEBUG = false; end;
    if isfield(ops, 'verbose'), VERBOSE = ops.verbose; else VERBOSE = true; end;
    if isfield(ops, 'f_fn'), f_fn = ops.f_fn; end;
    if isfield(ops, 'tol'), tol = ops.tol; else tol = 1e-6; end;
    if ~isfield(ops, 'prox_weight'), ops.prox_weight=1; end;

    if VERBOSE, fprintf('Running pgd_fista with DEBUG=%d ...\n', DEBUG); end;

    if ~isa(step_size,'function_handle') 
       stepsizeSchedule = @(i, g) step_size;
    else
       stepsizeSchedule = step_size;
    end

    tracking = struct;
    tracking.obj = 100*ones(max_iters, 1);
    tracking.norm_g2 = 100*ones(max_iters, 1);
    tracking.time = 100*ones(max_iters, 1);
    tracking.dis2 = 100*ones(max_iters, 1);
    tracking.Stepsize = 100*ones(max_iters, 1);
    tracking.nnz = 100*ones(max_iters, 1);

    duration = 0;
    step_tic = tic;
    x_km1 = init_point; y_k = x_km1;
    t_k = 1;
    duration = duration + toc(step_tic);
    for i=1:max_iters
        step_tic = tic;

        g = g_fn(y_k);
        stepSize = stepsizeSchedule(i, g);
        u = y_k - stepSize.*g;

        prox_weight = stepSize .* ops.prox_weight;
        x_k = p_fn(u, prox_weight);
        t_kp1 = (1 + sqrt(1+4*t_k^2))/2;
        y_kp1 = x_k + (t_k-1)/t_kp1*(x_k - x_km1);
        
        % stopping_criteria = norm(y_kp1 - y_k, 'fro');
        stopping_criteria = norm(x_k - x_km1, 'fro');
        x_km1 = x_k;
        y_k = y_kp1;
        t_k = t_kp1;

        duration = duration + toc(step_tic);
        if DEBUG, tracking.obj(i) = f_fn(x_k); end;
        if DEBUG, tracking.norm_g2(i) = trace(g'*g); end;
        if DEBUG, tracking.time(i) = duration; end;
        if DEBUG, dis_current = x_k - ops.ground_truth; end;
        if DEBUG, tracking.dis2(i) = trace(dis_current'*dis_current); end;
        if DEBUG, tracking.Stepsize(i) = mean(stepSize); end;
        if  stopping_criteria < tol, break; end;
        if ops.verbose
            if mod(i, 10) == 0
                fprintf(sprintf('Iter %d\n', i))
            end
        end


        tracking.nnz(i) = nnz(x_k);
    end
    tracking.num_iters = i;
    tracking.stopping_criteria = stopping_criteria;

    tracking.obj = tracking.obj(1:i);
    tracking.norm_g2 = tracking.norm_g2(1:i);
    tracking.time = tracking.time(1:i);
    tracking.dis2 = tracking.dis2(1:i);
    tracking.stopping_criteria = stopping_criteria;

end
