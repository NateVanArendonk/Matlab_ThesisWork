function pred_cor = quantile_correct( obs, obs_time, pred, pred_time )
%pred_cor = quantile_correct( obs, obs_time, pred, pred_time )
%   Corrects pred by applying quantile correction to CDF
%   Develops ECDFs for correction by using overlapping time series,
%   corection is applied to entire time series. Plotting show CDF's.


% Select just overlapping time series for estimating ECDFs
start_date = max([obs_time(1) pred_time(1)]);
end_date = min([obs_time(end) pred_time(end)]);

inds = obs_time < start_date | obs_time > end_date;
sync_obs_time = obs_time;
sync_obs_time(inds) = [];
sync_obs = obs;
sync_obs(inds) = [];

inds = pred_time < start_date | pred_time > end_date;
sync_pred_time = pred_time;
sync_pred_time(inds) = [];
sync_pred = pred;
sync_pred(inds) = [];

% Estimate ECDF for obs and pred, remove duplicated on x-axis
[fo,xo] = ecdf(sync_obs);
inds = find(diff(xo) == 0);
fo(inds) = [];
xo(inds) = [];
disp(length(inds))

[fp,xp] = ecdf(sync_pred);
inds = find(diff(xp) == 0);
fp(inds) = [];
xp(inds) = [];
disp(length(inds))

% Correct quantiles
% First find places on CDF where predicted values lie
q_probs = interp1(xp,fp,pred);
% Then use those CDF values to find where on obs cdf curve they should be,
pred_cor = interp1(fo,xo,q_probs);

% Quick plot
clf

% Plot original ecdfs
[fo,xo] = ecdf(obs);
[fp,xp] = ecdf(pred);
subplot(131)
hold on
plot(xo,fo,'.-')
plot(xp,fp,'.-')
grid on
legend('obs','pred')
title('Original time series')
ylabel('CDF')
xlabel('Variable')

% Plot overlapping time ecdfs
[fo,xo] = ecdf(sync_obs);
[fp,xp] = ecdf(sync_pred);
subplot(132)
hold on
plot(xo,fo,'.-')
plot(xp,fp,'.-')
grid on
title('Overlapping time series')
ylabel('CDF')
xlabel('Variable')

% Plot Corrected ecdfs
[fo,xo] = ecdf(obs);
[fp,xp] = ecdf(pred_cor);
subplot(133)
hold on
plot(xo,fo,'.-')
plot(xp,fp,'.-')
grid on
title('Quantile Corrected')
ylabel('CDF')
xlabel('Variable')


end

