##############################################
Maple code generating random Reuleaux heptagon
##############################################

restart;

sampleFloat := rand(0. .. 2); 
t1 := sampleFloat(); 
sampleFloat2 := rand(t1 .. 1.5); 
t2 := sampleFloat2(); 
sampleFloat3 := rand(t2 .. t2+.5); 
t3 := sampleFloat3(); 
sampleFloat4 := rand(t3 .. 3.1); 
t4 := sampleFloat4(); 
alpha := -sin(t2)+sin(t1)-sin(t4)+sin(t3); 
beta := 1+cos(t2)-cos(t1)+cos(t4)-cos(t3); 

t5 := -arcsin((1/2)*sqrt(alpha^2+beta^2))+arccos(alpha/sqrt(alpha^2+beta^2)); 
t6 := arcsin((1/2)*sqrt(alpha^2+beta^2))+arccos(alpha/sqrt(alpha^2+beta^2));

t1 := .2499109628; 
t2 := .7339481128; 
t3 := 1.205338126; 
t4 := 1.408149043; 
t5 := 1.875328710; 
t6 := 2.642468712;

##################################################
Constraint check to ensure closedness of the curve
##################################################

evalf(int(cos(t), t = t1 .. t2)+int(cos(t), t = t3 .. t4)+int(cos(t), t = t5 .. t6)); 
evalf(int(sin(t), t = t1 .. t2)+int(sin(t), t = t3 .. t4))+int(sin(t), t = t5 .. t6);

##############################################################################
Support function definition and plot of the Reuleaux heptagon and of its polar

                    7 pieces for h(t) over [0,pi]
##############################################################################

h1 := .5*cos(t); 
dh1 := diff(h1, t); 

h2 := h1+(int(cos(s), s = t1 .. t))*sin(t)-(int(sin(s), s = t1 .. t))*cos(t); 
dh2 := diff(h2, t); 

h3 := h1+(int(cos(s), s = t1 .. t2))*sin(t)-(int(sin(s), s = t1 .. t2))*cos(t); 
dh3 := diff(h3, t); 

h4 := h1+(int(cos(s), s = t1 .. t2)+int(cos(s), s = t3 .. t))*sin(t)-(int(sin(s), s = t1 .. t2)+int(sin(s), s = t3 .. t))*cos(t); 
dh4 := diff(h4, t); 

h5 := h1+(int(cos(s), s = t1 .. t2)+int(cos(s), s = t3 .. t4))*sin(t)-(int(sin(s), s = t1 .. t2)+int(sin(s), s = t3 .. t4))*cos(t); 
dh5 := diff(h5, t); 

h6 := h1+(int(cos(s), s = t1 .. t2)+int(cos(s), s = t3 .. t4)+int(cos(s), s = t5 .. t))*sin(t)-(int(sin(s), s = t1 .. t2)+int(sin(s), s = t3 .. t4)+int(sin(s), s = t5 .. t))*cos(t); 
dh6 := diff(h6, t); 

h7 := h1+(int(cos(s), s = t1 .. t2)+int(cos(s), s = t3 .. t4)+int(cos(s), s = t5 .. t6))*sin(t)-(int(sin(s), s = t1 .. t2)+int(sin(s), s = t3 .. t4)+int(sin(s), s = t5 .. t6))*cos(t); 
dh7 := diff(h7, t);

######### Plot of the heptagon (14 pieces) #########

g1 := plot([h1*cos(t)-dh1*sin(t), h1*sin(t)+dh1*cos(t), t = 0 .. t1]); 
g11 := plot([h1*cos(t)-dh1*sin(t)-cos(t), h1*sin(t)+dh1*cos(t)-sin(t), t = 0 .. t1]);  ### duplicate over [pi,2pi]
g2 := plot([h2*cos(t)-dh2*sin(t), h2*sin(t)+dh2*cos(t), t = t1 .. t2]); 
g21 := plot([h2*cos(t)-dh2*sin(t)-cos(t), h2*sin(t)+dh2*cos(t)-sin(t), t = t1 .. t2]); ### duplicate over [pi,2pi]
g3 := plot([h3*cos(t)-dh3*sin(t), h3*sin(t)+dh3*cos(t), t = t2 .. t3]); 
g31 := plot([h3*cos(t)-dh3*sin(t)-cos(t), h3*sin(t)+dh3*cos(t)-sin(t), t = t2 .. t3]); ### duplicate over [pi,2pi]
g4 := plot([h4*cos(t)-dh4*sin(t), h4*sin(t)+dh4*cos(t), t = t3 .. t4]); 
g41 := plot([h4*cos(t)-dh4*sin(t)-cos(t), h4*sin(t)+dh4*cos(t)-sin(t), t = t3 .. t4]); ### duplicate over [pi,2pi]
g5 := plot([h5*cos(t)-dh5*sin(t), h5*sin(t)+dh5*cos(t), t = t4 .. Pi]); 
g51 := plot([h5*cos(t)-dh5*sin(t)-cos(t), h5*sin(t)+dh5*cos(t)-sin(t), t = t4 .. t5]); ### duplicate over [pi,2pi]
g6 := plot([h6*cos(t)-dh6*sin(t), h6*sin(t)+dh6*cos(t), t = t5 .. t6]); 
g61 := plot([h6*cos(t)-dh6*sin(t)-cos(t), h6*sin(t)+dh6*cos(t)-sin(t), t = t5 .. t6]); ### duplicate over [pi,2pi]
g7 := plot([h7*cos(t)-dh7*sin(t), h7*sin(t)+dh7*cos(t), t = t6 .. Pi]); 
g71 := plot([h7*cos(t)-dh7*sin(t)-cos(t), h7*sin(t)+dh7*cos(t)-sin(t), t = t6 .. Pi]); ### duplicate over [pi,2pi]

plots[display]({g1, g11, g2, g21, g3, g31, g4, g41, g5, g51, g6, g61, g7, g71});

######### Plot of the polar set (14 pieces) #########

gp1 := plot([cos(t)/h1, sin(t)/h1, t = 0 .. t1]); 
gp11 := plot([cos(t)/(h1-1), sin(t)/(h1-1), t = 0 .. t1]); 
gp2 := plot([cos(t)/h2, sin(t)/h2, t = t1 .. t2]); 
gp21 := plot([cos(t)/(h2-1), sin(t)/(h2-1), t = t1 .. t2]); 
gp3 := plot([cos(t)/h3, sin(t)/h3, t = t2 .. t3]); 
gp31 := plot([cos(t)/(h3-1), sin(t)/(h3-1), t = t2 .. t3]); 
gp4 := plot([cos(t)/h4, sin(t)/h4, t = t3 .. t4]); 
gp41 := plot([cos(t)/(h4-1), sin(t)/(h4-1), t = t3 .. t4]); 
gp5 := plot([cos(t)/h5, sin(t)/h5, t = t4 .. t5]); 
gp51 := plot([cos(t)/(h5-1), sin(t)/(h5-1), t = t4 .. t5]); 
gp6 := plot([cos(t)/h6, sin(t)/h6, t = t5 .. t6]); 
gp61 := plot([cos(t)/(h6-1), sin(t)/(h6-1), t = t5 .. t6]); 
gp7 := plot([cos(t)/h7, sin(t)/h7, t = t6 .. Pi]); 
gp71 := plot([cos(t)/(h7-1), sin(t)/(h7-1), t = t6 .. Pi]);

plots[display]({gp1, gp11, gp2, gp21, gp3, gp31, gp4, gp41, gp5, gp51, gp6, gp61, gp7, gp71});


## Numerical attempt to compute the Santalo point and the Mahler area via maple (see also the Newton code in python)
fsolve(int(exp(I*t)/(h1-a*cos(t)-b*sin(t))^3, t = 0 .. t1)+int(exp(I*t)/(h2-a*cos(t)-b*sin(t))^3, t = t1 .. t2)+int(exp(I*t)/(h3-a*cos(t)-b*sin(t))^3, t = t2 .. t3)+int(exp(I*t)/(h4-a*cos(t)-b*sin(t))^3, t = t3 .. t4)+int(exp(I*t)/(h5-a*cos(t)-b*sin(t))^3, t = t4 .. t5)+int(exp(I*t)/(h6-a*cos(t)-b*sin(t))^3, t = t5 .. t6)+int(exp(I*t)/(h7-a*cos(t)-b*sin(t))^3, t = t6 .. t7)+int(exp(I*t)/(h7-a*cos(t)-b*sin(t))^3, t = t7 .. Pi) = 0, {a, b});
tmp1 := (1/2)*(int(2*h1*(h1+diff(h1, t$2)-1)+1, t = 0 .. t1))+(1/2)*(int(2*h2*(h2+diff(h2, t$2)-1)+1, t = 0*t1 .. t2))+(1/2)*(int(2*h3*(h3+diff(h3, t$2)-1)+1, t = t2 .. t3))+(1/2)*(int(2*h4*(h4+diff(h4, t$2)-1)+1, t = t3 .. t4))+(1/2)*(int(2*h5*(h5+diff(h5, t$2)-1)+1, t = t4 .. t5))+(1/2)*(int(2*h6*(h6+diff(h6, t$2)-1)+1, t = t5 .. t6))+(1/2)*(int(2*h7*(h7+diff(h7, t$2)-1)+1, t = t6 .. Pi))+subs(t = 0, diff(h1, t));



