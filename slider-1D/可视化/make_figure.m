  for i=1:6
    filetype='.dat';
    name=['E:\算例汇总\C-LCP-SQP-mSQP\slider-1D\可视化\sliderLCP',int2str(i)];
    filename=[name filetype];
    load(filename);
  end

figure(1)
plot(sliderLCP1(:,1),sliderLCP1(:,2),'b','linewidth',2);hold on
plot(sliderLCP1(:,1),sliderLCP1(:,3),'b-.','linewidth',2);hold on%SQP

figure(2)
plot(sliderLCP3(:,1),sliderLCP3(:,2),'k','linewidth',2);hold on
plot(sliderLCP3(:,1),sliderLCP3(:,3),'k-.','linewidth',2);hold on%Lemke

figure(3)
plot(sliderLCP5(:,1),sliderLCP5(:,2),'m','linewidth',2);hold on
plot(sliderLCP5(:,1),sliderLCP5(:,3),'m-.','linewidth',2);hold on%试算法

figure(4)
plot(sliderLCP2(:,1),sliderLCP2(:,2),'b','linewidth',2);hold on
plot(sliderLCP2(:,1),sliderLCP2(:,3),'b-.','linewidth',2);hold on%v_fsSQP
figure(5)
plot(sliderLCP4(:,1),sliderLCP4(:,2),'k','linewidth',2);hold on
plot(sliderLCP4(:,1),sliderLCP4(:,3),'k-.','linewidth',2);hold on%v_fsLCP
figure(6)
plot(sliderLCP6(:,1),sliderLCP6(:,2),'m','linewidth',2);hold on
plot(sliderLCP6(:,1),sliderLCP6(:,3),'m-.','linewidth',2);hold on%v_fs试算法
xlabel('时间t  [s]')
ylabel('相对速度、摩擦力  [m/s,N]')






