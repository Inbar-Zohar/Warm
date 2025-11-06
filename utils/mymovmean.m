function [t_movmean, v_movmean] = mymovmean(t,v,n)

v_movmean = movmean(v,n);
v_movmean = v_movmean(floor(n/2):n:end);
t_movmean = t(floor(n/2):n:end);