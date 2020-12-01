function r = weighted_corr(x,y,w)

x_mean = sum(sum(x.*w))./sum(sum(w));
y_mean = sum(sum(y.*w))./sum(sum(w));

cov_xy = sum(sum(w.*(x-x_mean).*(y-y_mean)))./sum(sum(w));
cov_xx = sum(sum(w.*(x-x_mean).*(x-x_mean)))./sum(sum(w));
cov_yy = sum(sum(w.*(y-y_mean).*(y-y_mean)))./sum(sum(w));

r=cov_xy/sqrt(cov_xx*cov_yy);

