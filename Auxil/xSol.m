function [x] = xSol(A, dSeq,t, xo, beta2 )
%	[x] = xSol(A, dSeq, t, xo, beta)
%		Computes the solution for the system with impulses given by dSeq. 
%		If xo is not given, the initial value is assumed to be zero.

	if (nargin < 5)
		beta2 = zeros(length(xo)-1,1);
	end
	B = [1;0];
	Beta = [0; beta2];
	[Td, Ad] = eig(A);
	Tdi = Td^(-1);
	n = length(A);

	x = zeros(n,length(t));
	%If xo is specified, there should be no impulses before t(1).
	if (nargin >= 4) 
		k = find(dSeq(1,:) <= t(1), 1, 'last');
		if isempty(k)
			k=0;
		end
		xi = xo;
		if (dSeq(1,1) == t(1))
			xi = xi+dSeq(2,1)*B;
		end
		dSeq = [[t(1); 0] dSeq(1:2, (k+1):end)];
	else
		dSeq = [ [min(dSeq(1,1)-10, t(1)-10);0] dSeq(1:2,:)];
		xi = zeros(n,1);
	end


	%Make sure that last there are firings both before and after t.
	if (dSeq(1,end) < t(end))
		dSeq = [dSeq [t(end)+10; 0]];
	end


	%xi is defined as x(t_i^+)

	for i=1:(length(dSeq)-1)
		from = find(t>=dSeq(1, i), 1, 'first');
		to = find(t<dSeq(1,i+1), 1, 'last');
		x(:, from:to) = expmtx(Ad,Td,Tdi, t(from:to)-dSeq(1,i), xi);
		xi = expmtx(Ad,Td,Tdi, dSeq(1,i+1) - dSeq(1,i), xi) + dSeq(2,i+1)*B;
	end

	x = x + A^-1*(expmtx(Ad, Td, Tdi, t,Beta) - expmtx(Ad,Td,Tdi,t*0,Beta));
