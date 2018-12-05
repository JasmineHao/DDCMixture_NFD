% Matrix Inversion

d = 1;
time_diff=0;
while time_diff < 1
     a = randn(d);
     ts = tic;
     inv(a);
     time_diff = toc(ts);
     time_used(d) = time_diff;
     d= d+1;
end

plot(1:length(time_used),time_used);
ylabel(['Time used(secons)']);
xlabel(['Matrix dimension'])
saveas(gcf,'figs/Matrix_Inversion.jpg');
