function [noise_fit] = fit_psd(f, psd)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% fit a PSD
fit_op = fitoptions('Method','NonlinearLeastSquares',...
                    'Lower',[0, 0, 0],...
                    'Upper',[Inf, Inf, Inf],...
                    'StartPoint',[2.013e-08, 2.715e-10, 5.135e-09]);

noise_model = fittype('(h2*x^2) + (h0) + (hm1/x)', ...
                      'coefficients', {'h2', 'h0', 'hm1',}, ...
                      'options', fit_op);

noise_fit = fit(f, psd, noise_model);


% % second try fit
fit_op = fitoptions('Method','NonlinearLeastSquares',...
                    'Lower',[0, 0, 0],...
                    'Upper',[Inf, Inf, Inf],...
                    'StartPoint',[noise_fit.h2, noise_fit.h0, ...
                    noise_fit.hm1]);
                
noise_model = fittype('(h2*x^2) + (h0) + (hm1/x)', ...
                      'coefficients', {'h2', 'h0', 'hm1',}, ...
                      'options', fit_op);

                
noise_fit = fit(f, psd, noise_model);



end



