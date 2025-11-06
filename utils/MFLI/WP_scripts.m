%% 
f=631;%631;%181;
DorB.Bx_sim.apply('Sin',f,0.1,0);
DorB.By_sim.apply('Sin',f,0.1,0);

%% Only By
DorB.Bx_sim.apply('OutputOFF')
DorB.By_sim.apply('OutputON')
%% Only Bx
DorB.Bx_sim.apply('OutputON')
DorB.By_sim.apply('OutputOFF')
%%
DorB.ESR_LIA.apply('shiftPhasePlus90')
%%
DorB.ESR_LIA.apply('shiftPhaseMinus90')
%%
DorB.ESR_LIA.apply('shiftPhase', -0.1)
%%
%% All Bsim on
DorB.Bx_sim.apply('OutputON')
DorB.By_sim.apply('OutputON')
%% All Bsim off
DorB.Bx_sim.apply('OutputOFF')
DorB.By_sim.apply('OutputOFF')


%%
%%
DorB.Xe129_LIA.apply('autoPhase')
%%
DorB.ESR_LIA.apply('shiftPhaseMinus90')

