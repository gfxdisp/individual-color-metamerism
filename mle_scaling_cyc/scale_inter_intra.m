function [phi_rec, delta_rec, v_rec, w_rec, OBSs, CONDs] = scale_inter_intra( observer_id, condition_id, quality, options )
arguments
    observer_id (:,1)
    condition_id (:,1)
    quality (:,1) { mustBeNumeric }
    options.no_delta (1,1) = false
end
% The function implememts scaling of MLE scaling of subjective judgements (qyuality) while accounting for
% within observer variance (v_rec) and inter-observer variance (w_rec).



if ~all(size(observer_id)==size(condition_id)) || ~all(size(condition_id)==size(quality))
    error( 'observer_id, condition_id and quality must be column vectors of the same size' );
end

[CONDs,~,cond_ind] = unique(condition_id);
N = numel(CONDs);
[OBSs,~,obs_ind] = unique(observer_id);
K = numel(OBSs);

phi_0 = ones(N,1) * mean(quality);
delta_0 = zeros(K,1);
log_v_0 = ones(K,1)*-2;
log_w_0 = ones(N,1)*-2;

if options.no_delta

    par_0 = [phi_0' log_v_0' log_w_0'];

    par_rec = fminunc( @log_likelyhood_no_delta, par_0 );

    phi_rec = par_rec(1:N);
    v_rec = exp(par_rec((N+1):(N+K)));
    w_rec = exp(par_rec((N+K+1):(2*N+K)));
    delta_rec = nan(size(v_rec));
else

    par_0 = [phi_0' delta_0' log_v_0'];
    par_rec = fminunc( @log_likelyhood, par_0 );

    phi_rec = par_rec(1:N);
    delta_rec = par_rec((N+1):(N+K));
    v_rec = exp(par_rec((N+K+1):(N+K*2)));

end

    function L = log_likelyhood( par )

        phi = par(1:N)';
        delta = par((N+1):(N+K))';
        log_v = par((N+K+1):(N+K*2))';
        log_w = par((N+K*2+1):(2*N+K*2))';

        log_vw = log(exp(log_v(obs_ind)) + exp(log_w(cond_ind)));

        L = sum( (quality-phi(cond_ind)-delta(obs_ind)).^2./(2*exp(log_vw).^2) + log_vw ); % + alpha*sum( delta.^2 );
    end

    function L = log_likelyhood_no_delta( par )

        phi = par(1:N)';
        log_v = par((N+1):(N+K))';
        log_w = par((N+K+1):(2*N+K))';

        log_vw = log(exp(log_v(obs_ind)) + exp(log_w(cond_ind)));

        L = sum( (quality-phi(cond_ind)).^2./(2*exp(log_vw).^2) + log_vw );
    end

end

