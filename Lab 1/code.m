
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3- Test the quantizer/dequantizer %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = -6:0.01:6;
n_bits= 3;
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
yline(0);
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
yline(0);
plot( x, deq_valcall2 );
%axis([-8 8 -8 8])
%axis padded
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc4_5 for req 4 and 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SNR_S_T = calc4_5(x, xmax, m)
SNR_S_T = zeros(2, 7);
n_bits = 2:1:8;
SNR_SIM = zeros(1, 7);
SNR_SIM_dB = zeros(1, 7);
SNR_THER = zeros(1, 7);
SNR_THER_dB = zeros(1, 7);
for i = 1:7
  q_indcall3 = UniformQuantizer(x, n_bits(i), xmax, m);
  deq_valcall3 = UniformDequantizer(q_indcall3, n_bits(i), xmax, m);
  error = deq_valcall3 - x;
  SNR_SIM(i) = mean(x.^2)/mean(error.^2);
  SNR_SIM_dB(i) = db(SNR_SIM(i));
  num_levels = 2 ^ n_bits(i);
  SNR_THER(i) = (mean(x.^2)*(3* num_levels^2 / xmax^2));
  SNR_THER_dB(i) = db(SNR_THER(i));

end
figure();
title("simulation and theoretical SNR DB")
hold on
plot(n_bits, SNR_SIM_dB,'r');
plot(n_bits, SNR_THER_dB,'b');
hold off
SNR_S_T(1,:) = SNR_SIM_dB;
SNR_S_T(2,:) = SNR_THER_dB;
% figure();
% title("simulation and theoretical SNR")
% hold on
% plot(n_bits, SNR_SIM,'r');
% plot(n_bits, SNR_THER,'b');
% hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  4- test your input on a random input signal  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = -5 + (5+5)*rand(10000,1);
xmax = 5;
m = 0;
calc4_5(x, xmax, m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  5- test your input on a random input signal  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_signal = exprnd(3,[1 10000]);
M=1;
N=10000;
random_sign_1 = mod( reshape(randperm(M*N), M, N), 2 );
random_sign = random_sign_1*2 - 1;
req_exp = random_sign.*exp_signal;
xmax = max(abs(req_exp));
m= 0;
SNR_S_T = calc4_5(req_exp, xmax, m);
signal_x = mean(req_exp.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%6- Now quantize the the non-uniform signal using a non-uniform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_signal = exprnd(3,[1 10000]);
M=1;
N=10000;
random_sign_1 = mod( reshape(randperm(M*N), M, N), 2 );
random_sign = random_sign_1*2 - 1;
req_exp = random_sign.*exp_signal;
xmax = max(abs(req_exp));
m= 0;
req_exp_normalized = req_exp / xmax;
mue = [5;100;200];
n_bits = 2:1:8;
SNR_SIM = zeros(1, 7);
SNR_SIM_dB = zeros(1, 7);
SNR_THER = zeros(1, 7);
SNR_THER_dB = zeros(1, 7);
colors_sim = ['b','g','m'];
colors_ther = ["--b","--g","--m"];
figure();
title("simulation and theoretical SNR DB")
hold on
plot(n_bits, SNR_S_T(1,:),'r');
plot(n_bits, SNR_S_T(2,:),'--r');
for k=1:numel(mue)
  compress = random_sign.*(log(1 + mue(k)*abs(req_exp_normalized))/log(1 + mue(k)));
  compress_max = max(abs(compress));
  for i = 1:7
      q_indcall6 = UniformQuantizer(compress, n_bits(i), compress_max, m);
      deq_valcall6 = UniformDequantizer(q_indcall6, n_bits(i), compress_max, m);
      expanded = random_sign .*(((1+mue(k)).^abs(deq_valcall6)-1)/mue(k));
      expanded_denormalized = expanded * xmax;
      error = expanded_denormalized - req_exp ;
      SNR_SIM(i) = mean(req_exp.^2)/mean(error.^2);
      SNR_SIM_dB(i) = db(SNR_SIM(i));
 
      num_levels = 2 ^ n_bits(i);
      SNR_THER(i) = 3*(num_levels)^2 / (log(1+mue(k)))^2 ;
      SNR_THER_dB(i) = db(SNR_THER(i));
  end
% figure();
% title("simulation and theoretical SNR DB","for mue = "+mue(k))
% hold on
plot(n_bits, SNR_SIM_dB,colors_sim(k));
plot(n_bits, SNR_THER_dB,colors_ther(k));
% hold off
% figure();
% title("simulation and theoretical SNR","for mue = "+mue(k))
% hold on
% plot(n_bits, SNR_SIM,colors_sim(k));
% plot(n_bits, SNR_THER,colors_ther(k));
% hold off
end
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1- UniformQuantizer Function%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q_ind = UniformQuantizer(in_val, n_bits, xmax, m)
  delta = deltaFunction(n_bits,xmax);
  q_ind = floor((in_val - ((m) * (delta / 2) - xmax)) / delta);
  q_ind(q_ind < 0) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2- UniformDequantizer Function%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deq_val = UniformDequantizer(q_ind, n_bits, xmax, m)
  delta = deltaFunction(n_bits,xmax);
  deq_val = ((q_ind) * delta) + ((m + 1) * (delta / 2) - xmax);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deltaFunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delta = deltaFunction(n_bits,xmax)
num_levels = 2 ^ n_bits; %2** n_bits;
delta = 2*xmax /num_levels;
end