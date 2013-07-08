function dy = myode(t, y, vt, v)

a1 = 97.8439;
a2 = -444.7561;
a3 = 112.2460;
a4 = -326.5010;

v = interp1(vt, v, t);

dy = zeros(6, 1);
dy(1) = y(2);
dy(2) = y(3);
dy(3) = y(4);

dy(4) = -(a1 + a2*v.*v).*y(2)-(a3 + a4*v.*v).*y;