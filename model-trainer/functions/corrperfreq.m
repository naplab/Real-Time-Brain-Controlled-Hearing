function [cr,cc] = corrperfreq (a,b)
% a(8:15,:)=[]; b(8:15,:)=[];
% a = a-mean2(a); b = b-mean2(b);
% a = a./std2(a);
% b = b./std2(b);
rng = 1:size(a,1);
% rng = 1:10;
for cnt10 = rng
    cc(cnt10) = corrnum(a(cnt10,:),b(cnt10,:));
%     cm(cnt10) = norm(a(cnt10,:)-b(cnt10,:));
end
% cr = mean(cm);
%% norm error:
% plot(cc);
cr = nanmean(cc);
% cr = norm(a-b);
% cr = mean(cc);
