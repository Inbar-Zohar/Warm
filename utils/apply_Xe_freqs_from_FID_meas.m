%%
DorB.Xe129Drive.apply('Sin',meas.FID_f_129,prm.WP.Xe129.amp,0)
DorB.Xe131Drive.apply('Sin',meas.FID_f_131,prm.WP.Xe131.amp,0)

return
%%
DorB.Xe129Drive.apply('Sin',prm.WP.Xe129.freq,prm.WP.Xe129.amp,0)
DorB.Xe131Drive.apply('Sin',prm.WP.Xe131.freq,prm.WP.Xe131.amp,0)
% 
% DorB.Xe129_PID.apply('setPropGain',-3)
% DorB.Xe129_PID.apply('setIntGain',0.3)
% 
% DorB.Xe131_PID.apply('setPropGain',-2)
% DorB.Xe131_PID.apply('setIntGain',0.08)
