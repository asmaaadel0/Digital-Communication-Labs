clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-Test the quantizer/dequantizer %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = -6:0.01:6;
n_bits= 2;
xmax = 6;
m = 1;
% midtread
q_indcall1 = UniformQuantizer(x, n_bits, xmax, m);
deq_valcall1 = UniformDequantizer(q_indcall1, n_bits, xmax, m);

figure();
title("Midtread")
hold on
plot( x, q_indcall1 );
plot( x, x );
plot( x, deq_valcall1 );
%axis([-8 8 -8 8])
hold off

m = 0;
% midrise
q_indcall2 = UniformQuantizer(x, n_bits, xmax, m);
deq_valcall2 = UniformDequantizer(q_indcall2, n_bits, xmax, m);

figure();
title("Midrise")
hold on
plot( x, q_indcall2 );
plot( x, x );
plot( x, deq_valcall2 );
%axis([-8 8 -8 8])
%axis padded
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UniformQuantizer Function%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function q_ind = UniformQuantizer(in_val, n_bits, xmax, m)
    num_levels = 2 ^ n_bits;
    delta = 2 * xmax / num_levels;
    q_ind = floor((in_val - ((m) * (delta / 2) - xmax)) / delta);
    q_ind(q_ind < 0) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UniformDequantizer Function%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function deq_val = UniformDequantizer(q_ind, n_bits, xmax, m)
    num_levels = 2 ^ n_bits;
    delta = 2 * xmax / num_levels;
    deq_val = ((q_ind) * delta) + ((m + 1) * (delta / 2) - xmax);
end

