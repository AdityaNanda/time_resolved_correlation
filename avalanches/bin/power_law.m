function [exp1, intrcpt,mse]=power_law(x,y)
% This function fits a power law of
% the form y= intrcpt.*x^exp1

x_guess= polyfit(log(x), log(y),1); 
exp1= x_guess(1); intrcpt= exp(x_guess(2)); 


 w=ones(size(x)); %.* 1./x.^x_guess(1); 
 func= @(var) w.*(y- var(2).*x.^var(1)); 
 mse= norm(func((x_guess))).^2; 
% opt=optimoptions('lsqnonlin', 'display', 'iter', ...
%     'algorithm', 'levenberg-marquardt',...
%     'MaxIter',1e3, 'FunctionTolerance', 1e-16 ); 
% [x0, resnorm, res] = lsqnonlin(func, (x_guess),[],[], opt); 
% w=ones(size(x)); 
% res= y- x0(1).*x.^x0(2); 
% mse= sqrt(resnorm);
% exp1= x0(2); 
% intrcpt= x0(1); 

end


% close all; 
% loglog(x,y, 'r', 'linewidth', 2); 
% hold on; 
% loglog(x,intrcpt.*x.^exp1, 'linewidth',2)
% xlabel('T'); 
% ylabel('Average size (T), <s>(T)' ); 
% set(gca, 'fontsize',27); 
% legend('values', 'nonlinlsq fit'); 
% print([plot_save_folder, 'plots/fit_lin'], '-dpng', '-r340'); 