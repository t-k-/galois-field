% example 4.6 at page 97
% We want to multiply the two polynomials:
% A(x) = x^3 + x^2 + 1 and B(x) = x^2 + x in the field GF(2^4)
% the irreducible polynomial of this Galois field is given as:
% P(x) = x^4 + x + 1
% The plain polynomial product is computed as:
% C(x) = A(x) · B(x) = x^5 + x^3 + x^2 + x
% C(x) = x^3 mod P(x)

%Galois_all_irreducibles(2, 4)
%Galois_test_irreducible([1 0 0 0 1 1], 2)
%Galois_test_irreducible([1 0 0 1 0 1], 2)
%Galois_test_irreducible([1 0 1 0 1 1], 2)
%Galois_all_irreducibles(3, 3)

%Ax = [1 0 1];
%Bx = [1 0 1 1];
%Ax = [1 0 1];
%Bx = [1 1];
%prim_poly = [1 1 0 0 1];
%polystring(Galois_rawmul(Ax, Bx, 2), '\n')
%Galois_visualmul(Ax, Bx, 2, prim_poly)

% Ax = [1 1 0 1 1];
% Bx = [0 0 1 1 0 0 1 1];
% prim_poly = [1 1 0 1 1 0 0 0 1];
% Galois_visualmul(Ax, Bx, 2, prim_poly)

%Ax = [1 1 0 0 1];
%Bx = [1 1 0 1 1];
%prim_poly = [1 1 0 1 1 0 0 0 1];
%polystring(Galois_rawmul(Ax, Bx, 2), '\n')
%Galois_visualmul(Ax, Bx, 2, prim_poly)

% Ax = [1 1 0 0 1 0 1 1];
% Bx = [1 1 0 0 0 1 0 1];
% prim_poly = [1 1 0 1 1 0 0 0 1];
% polystring(Galois_rawmul(Ax, Bx, 2), '\n')
% Galois_visualmul(Ax, Bx, 2, prim_poly)

prim_poly = [1 1 0 0 1];
fprintf('~~~ generating using polynomial: '); polystring(prim_poly, ' ~~~\n')
Galois_pow([0 1], 23, 2, prim_poly)

%   Ax = [1 0 1 1];
%   Bx = [0 1 1];
%   prim_poly = [1 1 0 0 1];
%   Galois_visualmul(Ax, Bx, 2, prim_poly)
%   
%   % non-monic casting to monic
%   Ax = [2 2];
%   Bx = [2];
%   prim_poly = [1 0 1];
%   Galois_visualmul(Ax, Bx, 3, prim_poly)
%   
%   Ax = [1 0 0 1];
%   Bx = [1 1 1];
%   Galois_visualrawdiv(Ax, Bx, 2)
%   
%   Ax = [1 0 1 1];
%   Bx = [1 1];
%   Galois_visualrawdiv(Ax, Bx, 2)
%   
%   Ax = [1 2 1];
%   Bx = [1 1];
%   Galois_visualrawdiv(Ax, Bx, 3)
%   
%   p = 2;
%   for m = 0:4
%   	fprintf('~~~ base = %d, degree = %d ~~~\n', p, m)
%   	Galois_all_irreducibles(p, m)
%   end
%   

% p = 3;
% for m = 0:2
% 	fprintf('~~~ base = %d, degree = %d ~~~\n', p, m)
% 	Galois_all_irreducibles(p, m)
% end
% 
% prim_poly = [2 2 1];
% fprintf('~~~ generating using polynomial: '); polystring(prim_poly, ' ~~~\n')
% Galois_generate([1 1], 3, 2, prim_poly)
% 
% prim_poly = [2 1 1];
% fprintf('~~~ generating using polynomial: '); polystring(prim_poly, ' ~~~\n')
% Galois_generate([1 1], 3, 2, prim_poly)
% 
% prim_poly = [1 0 1];
% fprintf('~~~ generating using polynomial: '); polystring(prim_poly, ' ~~~\n')
% Galois_generate([1 1], 3, 2, prim_poly)

%Galois_visualrawdiv([1 1 0 0 1], [0 0 1 1 0 0 1 1], 2)

function Galois_pow(Ax, power, p, prim_poly)
	Cx = [1];
	for pow = 0: power
		fprintf('(')
		polystring(Ax, ')^')
		fprintf('%d', pow)
		fprintf(' = ')
		polystring(Cx, '\n')
		Cx = Galois_mul(Ax, Cx, p, prim_poly);
	end
end

function Galois_generate(Ax, p, m, prim_poly)
	Cx = [1];
	seen = cell(0);
	n_elem = p^m;
	if n_elem < 2
		return
	end
	for pow = 0: n_elem - 2
		fprintf('(')
		polystring(Ax, ')^')
		fprintf('%d', pow)
		fprintf(' = ')
		polystring(Cx, '\n')

		found = false;
		for idx = 1:length(seen)
			if seen{idx} == Cx
				if length(seen{idx}) == length(Cx)
					found = true;
					disp('start to repeat ...')
					break;
				end
			end
		end

		if found
			break;
		end
		seen = [seen; Cx];
		Cx = Galois_mul(Ax, Cx, p, prim_poly);
	end
end

function Galois_test_irreducible(Ax, base)
	degree = find(Ax, 1, 'last') - 1;
	[candidates, denominators] = Galois_irr_candidates(base, degree);
	[j_rows, discard_] = size(denominators);

	is_irreducible = true;
	for j = 1:j_rows
		j_row = denominators(j, :);
		j_row = j_row(1:find(j_row, 1, 'last'));
		Galois_visualrawdiv(Ax, j_row, base);
		[Qx Rx] = Galois_rawdiv(Ax, j_row, base);
		if length(Rx) == 0
			is_irreducible = false;
			break
		end
	end
	if is_irreducible
		polystring(Ax, ' ');
		disp('is irreducible !!!');
	end
end

function Galois_all_irreducibles(base, degree)
	[candidates, denominators] = Galois_irr_candidates(base, degree);
	[i_rows, discard_] = size(candidates);
	[j_rows, discard_] = size(denominators);
	for i = 1:i_rows
		i_row = candidates(i, :);
		if degree > 0 && i_row(1) == 0
			continue
		end
		is_irreducible = true;
		for j = 1:j_rows
			j_row = denominators(j, :);
			j_row = j_row(1:find(j_row, 1, 'last'));
			Galois_visualrawdiv(i_row, j_row, base);
			[Qx Rx] = Galois_rawdiv(i_row, j_row, base);
			if length(Rx) == 0
				is_irreducible = false;
				break
			end
		end
		if is_irreducible
			polystring(i_row, ' ');
			disp('is irreducible !!!');
		end
	end
end

function [candidates, denominators] = Galois_irr_candidates(base, degree)
	all_poly = perm_comb(base, degree + 1);
	[rows, columns] = size(all_poly);
	candidate_prim_polys = [];
	denominator_polys = [];
	for i = 1:rows
		row = all_poly(i, :);
		if row == zeros(1, degree + 1)
			if degree == 0
				candidate_prim_polys = [row; candidate_prim_polys];
			end
			continue
		end
		row_deg = find(row, 1, 'last') - 1;
		rightmost_coeff = row(row_deg + 1);
		if row_deg == degree
			% filter out non-monic
			if rightmost_coeff == 1 || degree == 0
				candidate_prim_polys = [row; candidate_prim_polys];
			end
		elseif row_deg ~= 0
			% filter out non-monic
			if rightmost_coeff == 1
				denominator_polys = [row; denominator_polys];
			end
		end
	end
	candidates = candidate_prim_polys;
	denominators = denominator_polys;
end

function Galois_visualmul(Ax, Bx, p, prim_poly)
	fprintf('(')
	polystring(Ax, ') * (')
	polystring(Bx, ') = ')
	Cx = Galois_mul(Ax, Bx, p, prim_poly);
	polystring(Cx, '  (using primitive poly: ')
	polystring(prim_poly, ')\n')
end

function Galois_visualrawdiv(Ax, Bx, p)
	polystring(Ax, ' = (')
	polystring(Bx, ') * (')
	[Qx Rx] = Galois_rawdiv(Ax, Bx, p);
	polystring(Qx, ') + (')
	polystring(Rx, ') \n')
end

function CPrime = Galois_rawmul(Ax, Bx, p)
	polynomials_mul = conv(fliplr(Ax), fliplr(Bx)); % polynomials multiplication
	polynomials_mod = mod(polynomials_mul, p); % modulo p
	CPrime = fliplr(polynomials_mod); % flip vector
	CPrime = CPrime(1:find(CPrime, 1, 'last')); % cut trailing zeros
end

function Cx = Galois_mul(Ax, Bx, p, prim_poly)
	CPrime = Galois_rawmul(Ax, Bx, p);
	[q r] = deconv(fliplr(CPrime), fliplr(prim_poly));
	Cx = fliplr(mod(r, p));
	Cx = Cx(1:find(Cx, 1, 'last'));
end

function [Qx, Rx] = Galois_rawdiv(Ax, Bx, p)
	[q r] = deconv(fliplr(Ax), fliplr(Bx));
	Qx = fliplr(mod(q, p));
	Qx = Qx(1:find(Qx, 1, 'last')); % cut trailing zeros
	Rx = fliplr(mod(r, p));
	Rx = Rx(1:find(Rx, 1, 'last')); % cut trailing zeros
end

function polystring(Ax, endstr)
	stack = [];
	for i = 1:length(Ax)
		ord_str  = num2str(i - 1);
		coef_str = num2str(Ax(i));
		if Ax(i) ~= 0
			if i == 1
				stack = [stack; coef_str];
			else
				if Ax(i) == 1
					coef_str = '';
				end
				if i - 1 == 1
					stack = [stack; string(strcat(coef_str, 'x'))];
				else
					stack = [stack; string(strcat(coef_str, 'x^', ord_str))];
				end
			end
		end
	end
	stack = flipud(stack);
	stacklen = length(stack);
	if stacklen == 0
		fprintf('0')
	end
	for i = 1:stacklen
		if i == stacklen
			fprintf('%s', stack(i));
		else
			fprintf('%s + ', stack(i));
		end
	end
	fprintf(endstr)
end

function A = perm_comb(p, m)
	A = combinator(p, m, 'p', 'r');
	A(A == p) = 0;
end
