function Zeta()

    a1 = 97.8439;
    a2 = -444.7561;
    a3 = 112.2460;
    a4 = -326.5010;
    a1s = 0.4519;
    a2s = 0.3317;

    v_o = 1000;
    v_D = 200;
    OmegaBar = 2;
    
    v = @(t) (v_o + v_D*cos(OmegaBar*t));
    
    [T, Y] = ode45(@linearized, [0 200], [sqrt(2) 0 0 0]);
    
    
    xlabel('\tau')
    ylabel('\zeta');
    plot(T, Y(:,1), '-k')
    
    function dy = linearized(t,y)
       dy = zeros(4,1);
       dy(1) = y(2);
       dy(2) = y(3);
       dy(3) = y(4);
       %dy(4) = - (a1+random('Normal',0, 1)*random('Normal',0, 1)*random('Normal',0, 1))*y(3) - (a3+random('Normal',0, 1)*random('Normal',0, 1)*random('Normal',0, 1))*y(1);
       dy(4) = -a1s*v(t)*y(4) - (a1)*y(3) - a2s*v(t)*y(2) - (a3)*y(1);
    end    
end
