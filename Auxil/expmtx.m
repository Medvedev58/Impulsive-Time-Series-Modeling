function [eAT] = expmtx(Ad, Td, Tdi, t,x0)
% expmtx(Ad,Td,Tdi,t,x0)	Computes expm(A*t)*x0 for all times in t, where
%				A = Td*Ad*Tdi, Ad is diagonal and Tdi = Td^-1
    w0 = Tdi*x0;

    w = [];

	
    for i = 1:length(x0)
    	wi = exp(Ad(i,i)*t)*w0(i);
		w = [w; wi];
    end

    eAT = Td*w;

end
