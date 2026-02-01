function [c,s] = corrnum(x,y,flag)
if ~exist('flag','var') || isempty(flag)
    flag = 0;
end
s = []; c = [];
x=x(:);y=y(:);
if ~flag
    xm=mean(x);ym=mean(y);
    c=(mean(x.*y)-xm*ym)/sqrt(mean((x-xm).^2)*mean((y-ym).^2));
else
    %20?
    N = 10;
    sz = floor(size(x,1)/N);
%     ind = randperm(length(x));
%     x = x(ind);y=y(ind);
    for cnt2 = 1:N
        tst = (cnt2-1)*sz+1 : cnt2*sz;
        if cnt2 == N, tst(end+1:end+size(x,1)-tst(end)) = tst(end)+1:size(x,1);end
        xt = x(tst); yt = y(tst);
        xm=mean(xt);ym=mean(yt);
        c(cnt2) = (mean(xt.*yt)-xm*ym)/sqrt(mean((xt-xm).^2)*mean((yt-ym).^2));
    end
    s = std(c)/sqrt(N);
    c = mean(c);
end


