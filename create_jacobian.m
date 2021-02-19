function js = create_jacobian()
    syms H12 H12y1 H12y2
    syms H34 H34y1 H34y2
    syms H36 H36y1 H36y2
    syms H52 H52y1 H52y2
    syms H23 H23s1 H23s2 H23y1 H23y2
    syms H41 H41s1 H41s2 H41y1 H41y2
    syms H65 H65s1 H65s2 H65y1 H65y2

    syms IhvAw IlvAw IhvBw IlvBw IhvCw IlvCw
    syms delUhvAw delUlvAw delUhvBw delUlvBw delUhvCw delUlvCw
    syms UhvA UhvB UhvC UlvA UlvB UlvC
    syms IhvA IhvB IhvC IlvA IlvB IlvC
    syms RhvA RhvB RhvC RlvA RlvB RlvC
    syms Un Rn

    syms l12 l12y1 l12y2
    syms l23 l23y1 l23y2 l23s1 l23s2
    syms l34 l34y1 l34y2
    syms l36 l36y1 l36y2
    syms l41 l41y1 l41y2 l41s1 l41s2
    syms l52 l52y1 l52y2
    syms l65 l65y1 l65y2 l65s1 l65s2

    syms s12 s23 s34 s36 s41 s52 s65

    syms WhvA WlvA WhvB WlvB WhvC WlvC

    mu0 = 4 * pi * 10^(-7);

    f_ferr = @(H, s)(mu0 * sigmoid(H, 1, 1) * H * s);
    f_air = @(H, s)(mu0 * H * s);

    H_EXPR1  = H41s1 * l41s1 - H41 * l41 - IhvAw * WhvA + IlvAw * WlvA;
    H_EXPR2  = H41y1 * l41y1 - H41 * l41 - IhvAw * WhvA + IlvAw * WlvA;
    H_EXPR3  = H41y2 * l41y2 - H41 * l41 - IhvAw * WhvA + IlvAw * WlvA;
    H_EXPR4  = H41s2 * l41s2 - H41 * l41 - IhvAw * WhvA + IlvAw * WlvA;
    H_EXPR5  = H23s1 * l23s1 - H23 * l23 + IhvBw * WhvB - IlvBw * WlvB;
    H_EXPR6  = H23y1 * l23y1 - H23 * l23 + IhvBw * WhvB - IlvBw * WlvB;
    H_EXPR7  = H23y2 * l23y2 - H23 * l23 + IhvBw * WhvB - IlvBw * WlvB;
    H_EXPR8  = H23s2 * l23s2 - H23 * l23 + IhvBw * WhvB - IlvBw * WlvB;
    H_EXPR9  = H65s1 * l65s1 - H65 * l65 - IhvCw * WhvC + IlvCw * WlvC;
    H_EXPR10 = H65y1 * l65y1 - H65 * l65 - IhvCw * WhvC + IlvCw * WlvC;
    H_EXPR11 = H65y2 * l65y2 - H65 * l65 - IhvCw * WhvC + IlvCw * WlvC;
    H_EXPR12 = H65s2 * l65s2 - H65 * l65 - IhvCw * WhvC + IlvCw * WlvC;
    H_EXPR13 = H41 * l41 + H12 * l12 + H23 * l23 + H34 * l34 ...
        - IhvAw * WhvA + IhvBw * WhvB + IlvAw * WlvA - IlvBw * WlvB;
    H_EXPR14 = H23 * l23 + H36 * l36 + H65 * l65 + H52 * l52 ...
        + IhvBw * WhvB - IhvCw * WhvC - IlvBw * WlvB + IlvCw * WlvC;
    H_EXPR15 = H12y2 * l12y2 - H12 * l12;
    H_EXPR16 = H12y1 * l12y1 - H12 * l12;
    H_EXPR17 = H52y2 * l52y2 - H52 * l52;
    H_EXPR18 = H52y1 * l52y1 - H52 * l52;
    H_EXPR19 = H34y2 * l34y2 - H34 * l34;
    H_EXPR20 = H34y1 * l34y1 - H34 * l34;
    H_EXPR21 = H36y2 * l36y2 - H36 * l36;
    H_EXPR22 = H36y1 * l36y1 - H36 * l36;

    F_EXPR1 = f_ferr(H41, s41) + f_air(H41s1, s41) + f_air(H41s2, s41) ...
        + f_air(H41y1, s41) + f_air(H41y2, s41) + f_ferr(H12, s12) ...
        + f_air(H12y1, s12) - f_air(H12y2, s12);
    F_EXPR2 = f_ferr(H52, s52) + f_air(H52y1, s52) + f_air(H52y2, s52) ...
        + f_ferr(H12, s12) + f_air(H12y1, s12) + f_air(H12y2, s12) ...
        - f_ferr(H23, s23) - f_air(H23s1, s23) - f_air(H23s2, s23) ...
        - f_air(H23y1, s23) - f_air(H23y2, s23);
    F_EXPR3 = f_ferr(H65, s65) + f_air(H65s1, s65) + f_air(H65s2, s65) ...
        + f_air(H65y1, s65) + f_air(H65y2, s65) - f_ferr(H52, s52) ...
        - f_air(H52y1, s52) - f_air(H52y2, s52);
    F_EXPR4 = f_ferr(H36, s36) + f_air(H36y1, s36) + f_air(H36y2, s36) ...
        - f_ferr(H65, s65) - f_air(H65s1, s65) - f_air(H65s2, s65) ...
        - f_air(H65y1, s65) - f_air(H65y2, s65);
    F_EXPR5 = f_ferr(H23, s23) + f_air(H23s1, s23) + f_air(H23s2, s23) ...
        + f_air(H23y1, s23) + f_air(H23y2, s23) - f_ferr(H34, s34) ...
        - f_air(H34y1, s34) - f_air(H34y2, s34) - f_ferr(H36, s36) ...
        - f_air(H36y1, s36) - f_air(H36y2, s36);

    U_EXPR1 = delUhvAw - IhvAw * RhvA;
    U_EXPR2 = delUhvBw - IhvBw * RhvB;
    U_EXPR3 = delUhvCw - IhvCw * RhvC;
    U_EXPR4 = delUlvAw - IlvAw * RlvA;
    U_EXPR5 = delUlvBw - IlvBw * RlvB;
    U_EXPR6 = delUlvCw - IlvCw * RlvC;

    U_EXPR7  = delUhvAw - UhvA + Un;
    U_EXPR8  = delUhvBw - UhvB + Un;
    U_EXPR9  = delUhvCw - UhvC + Un;
    U_EXPR10 = Un - (IhvAw + IhvBw + IhvCw) * Rn;

    U_EXPR11 = delUlvAw - UlvA + UlvC;
    U_EXPR12 = delUlvBw - UlvB + UlvA;
    U_EXPR13 = delUlvCw - UlvC + UlvB;

    I_EXPR1 = IhvA - IhvAw;
    I_EXPR2 = IhvB - IhvBw;
    I_EXPR3 = IhvC - IhvCw;

    I_EXPR4 = IlvA - IlvAw + IlvBw;
    I_EXPR5 = IlvB - IlvBw + IlvCw;
    I_EXPR6 = IlvC - IlvCw + IlvAw;

    j = jacobian([H_EXPR1, H_EXPR2, H_EXPR3, H_EXPR4, H_EXPR5, H_EXPR6, ...
        H_EXPR7, H_EXPR8, H_EXPR9, H_EXPR10, H_EXPR11, H_EXPR12, ...
        H_EXPR13, H_EXPR14, H_EXPR15, H_EXPR16, H_EXPR17, H_EXPR18, ...
        H_EXPR19, H_EXPR20, H_EXPR21, H_EXPR22, F_EXPR1, F_EXPR2, ...
        F_EXPR3, F_EXPR4, F_EXPR5, U_EXPR1, U_EXPR2, U_EXPR3, U_EXPR4, ...
        U_EXPR5, U_EXPR6, U_EXPR7, U_EXPR8, U_EXPR9, U_EXPR10, ...
        U_EXPR11, U_EXPR12, U_EXPR13, I_EXPR1, I_EXPR2, I_EXPR3, ...
        I_EXPR4, I_EXPR5, I_EXPR6], ...
        [H12, H12y1, H12y2, H34, H34y1, H34y2, H36, H36y1, H36y2, H52, ...
        H52y1, H52y2, H23, H23s1, H23s2, H23y1, H23y2, H41, H41s1, ...
        H41s2, H41y1, H41y2, H65, H65s1, H65s2, H65y1, H65y2, IhvA, ...
        IhvB, IhvC, IlvA, IlvB, IlvC, RhvA, RhvB, RhvC, RlvA, RlvB, ...
        RlvC, Rn]);
    
    js = subs(j, ...
        {l12, l12y1, l12y2, l23, l23y1, l23y2, l23s1, l23s2, l34, ...
        l34y1, l34y2, l36, l36y1, l36y2, l41, l41y1, l41y2, l41s1, ...
        l41s2, l52, l52y1, l52y2, l65, l65y1, l65y2, l65s1, l65s2, ...
        s12, s23, s34, s36, s41, s52, s65, WhvA, WlvA, WhvB, WlvB, ...
        WhvC, WlvC}, ...
        {0.5, 0.5 * 1.2, 0.5 * 1.2, 0.8, 0.8 * 1.2, 0.8 * 1.2, ...
        0.8 * 1.3, 0.8 * 1.3, 0.5, 0.5 * 1.2, 0.5 * 1.2, 0.5, ...
        0.5 * 1.2, 0.5 * 1.2, 0.8, 0.8 * 1.2, 0.8 * 1.2, 0.8 * 1.3, ...
        0.8 * 1.3, 0.5, 0.5 * 1.2, 0.5 * 1.2, 0.8, 0.8 * 1.2, ...
        0.8 * 1.2, 0.8 * 1.3, 0.8 * 1.3, 64 * 10^(-4), 64 * 10^(-4), ...
        64 * 10^(-4), 64 * 10^(-4), 64 * 10^(-4), 64 * 10^(-4), ...
        64 * 10^(-4), 350, 100, 350, 100, 350, 100});
end