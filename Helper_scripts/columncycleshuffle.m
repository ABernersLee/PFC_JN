function [wc,p] = columncycleshuffle(post,numshuff)

extend_x = repmat(1:size(post,1),size(post,2),1);
extend_y = repmat([1:size(post,2)]',1,size(post,1));
wc = weighted_corr(extend_x,extend_y,post'); 
pp = randi(size(post,1),size(post,1),numshuff);

wcs = NaN(numshuff,1);
% tic
for i = 1:numshuff
    extend_x = repmat(1:size(post(pp(:,i),:),1),size(post(pp(:,i),:),2),1);
    extend_y = repmat([1:size(post(pp(:,i),:),2)]',1,size(post(pp(:,i),:),1));
    wcs(i) = weighted_corr(extend_x,extend_y,post(pp(:,i),:)'); 
end
% toc

p = (sum(abs(wcs)>=abs(wc))+1)/(numshuff+1);