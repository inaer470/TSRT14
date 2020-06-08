function [y] = y_TDOA2(tphat)
for m = 1:length(tphat)
     for i =[1:6],
        yy(m,i) = tphat(m,7)-tphat(m,i);
     end
end

y = sig(yy, 2);
end

