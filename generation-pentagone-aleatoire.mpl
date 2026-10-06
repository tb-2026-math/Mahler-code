

restart;
sampleFloat := rand(0. .. 2); t1 := sampleFloat(); sampleFloat2 := rand(t1 .. 3.1); t2 := sampleFloat2(); alpha := -sin(t2)+sin(t1); beta := 1+cos(t2)-cos(t1); t3 := -arcsin((1/2)*sqrt(alpha^2+beta^2))+arccos(alpha/sqrt(alpha^2+beta^2)); t4 := arcsin((1/2)*sqrt(alpha^2+beta^2))+arccos(alpha/sqrt(alpha^2+beta^2));
evalf(int(cos(t), t = t1 .. t2)+int(cos(t), t = t3 .. t4)); evalf(int(sin(t), t = t1 .. t2)+int(sin(t), t = t3 .. t4));
h1 := .5*cos(t); dh1 := diff(h1, t); h2 := h1+(int(cos(s), s = t1 .. t))*sin(t)-(int(sin(s), s = t1 .. t))*cos(t); dh2 := diff(h2, t); h3 := h1+(int(cos(s), s = t1 .. t2))*sin(t)-(int(sin(s), s = t1 .. t2))*cos(t); dh3 := diff(h3, t); h4 := h1+(int(cos(s), s = t1 .. t2)+int(cos(s), s = t3 .. t))*sin(t)-(int(sin(s), s = t1 .. t2)+int(sin(s), s = t3 .. t))*cos(t); dh4 := diff(h4, t); h5 := h1+(int(cos(s), s = t1 .. t2)+int(cos(s), s = t3 .. t4))*sin(t)-(int(sin(s), s = t1 .. t2)+int(sin(s), s = t3 .. t4))*cos(t); dh5 := diff(h5, t);

