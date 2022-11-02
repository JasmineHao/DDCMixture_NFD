%%

x = optimvar('x');
y = optimvar('y');

obj = objfunx(x,y);
prob = optimproblem('Objective',obj);

x0.x = -1;
x0.y = 1;
[sol2,fval2] = solve(prob,x0)

%%
TiltEllipse = x.*y/2 + (x+2).^2 + (y-2).^2/2 <= 2;
prob.Constraints.constr = TiltEllipse;
x0.x = -3;
x0.y = 3;
show(prob)
[sol,fval] = solve(prob,x0)

%%
x0.x = -1;
x0.y = 1;
[sol2,fval2] = solve(prob,x0)
%%
f = @objfunx;
g = @(x,y) x.*y/2+(x+2).^2+(y-2).^2/2-2;
rnge = [-5.5 -0.25 -0.25 7];
fimplicit(g,'k-')
axis(rnge);
hold on
fcontour(f,rnge,'LevelList',logspace(-1,1))
plot(sol.x,sol.y,'ro','LineWidth',2)
plot(sol2.x,sol2.y,'ko','LineWidth',2)
legend('Constraint','f Contours','Global Solution','Local Solution','Location','northeast');
hold off