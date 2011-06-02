function plot_Lhat_v_n(keeps,ns,fname,savestuff)

figure(1), clf

gray1=0.5*[1 1 1];
gray2=0.75*[1 1 1];
sa=0.5;
fs=12;


%
h(1)=subplot(211); hold all
errorbar(ns,keeps(1).misedge_rate,keeps(1).misedge_ste,'k','linewidth',2)
errorbar(ns,keeps(2).misedge_rate,keeps(2).misedge_ste,'color',gray1,'linewidth',2)
axis([0 ns(end) 0 1])
% set(gca,'YTick',[2:2:s])
xlabel('# training samples','fontsize',fs)
ylabel('missed-edge rate','fontsize',fs)

h(2)=subplot(212); hold all
errorbar(ns,keeps(1).Lhat_avg,keeps(1).Lhat_ste,'k','linewidth',2)
errorbar(ns,keeps(2).Lhat_avg,keeps(2).Lhat_ste,'color',gray1,'linewidth',2)
errorbar(ns,keeps(3).Lhat_avg,keeps(3).Lhat_ste,'color',gray2,'linewidth',2)
axis([0 ns(end) 0 0.5])
xlabel('# training samples','fontsize',fs)
ylabel('misclassification rate','fontsize',fs)
% legend(alg(1).name,alg(2).name, alg(3).name,'fontsize',fs)
legend('coh','inc','nb'); %,'fontsize',fs,'Location','Best')
set(gca,'YTick',[0.1:0.1:0.5])



% for i=1:numel(h)
%     op=get(h(i),'OuterPosition');
%     op(4)=op(4)+op(2);
%     op(2)=0;
%     set(h(i),'OuterPosition',op)
%     set(h(i),'fontsize',fs)
% end

if savestuff==1
    print_fig(['../../figs/' fname '_Lhats'],[3 3]*1.5)
end

