  for i=1:6
    filetype='.dat';
    name=['C:\Users\Administrator\Desktop\test\sliderLCP',int2str(i)];
    filename=[name filetype];
    load(filename);
  end


plot(sliderLCP1(:,1),sliderLCP1(:,2),'b');hold on
plot(sliderLCP1(:,1),sliderLCP1(:,3),'b-.');hold on%SQP

plot(sliderLCP3(:,1),sliderLCP3(:,2),'k','linewidth',2);hold on
plot(sliderLCP3(:,1),sliderLCP3(:,3),'k-.','linewidth',2);hold on%LCP

plot(sliderLCP5(:,1),sliderLCP5(:,2),'m','linewidth',3);hold on
plot(sliderLCP5(:,1),sliderLCP5(:,3),'m-.','linewidth',3);hold on%LCP_ref

plot(sliderLCP2(:,1),sliderLCP2(:,2),'b','linewidth',2);hold on%v_fs
plot(sliderLCP2(:,1),sliderLCP2(:,3),'b-.','linewidth',2);hold on%v_fsSQP

plot(sliderLCP4(:,1),sliderLCP4(:,2),'k','linewidth',2);hold on%v_fs
plot(sliderLCP4(:,1),sliderLCP4(:,3),'k-.','linewidth',2);hold on%v_fsLCP

plot(sliderLCP6(:,1),sliderLCP6(:,2),'m','linewidth',3);hold on%v_fs
plot(sliderLCP6(:,1),sliderLCP6(:,3),'m-.','linewidth',3);hold on%v_fsLCP-ref






