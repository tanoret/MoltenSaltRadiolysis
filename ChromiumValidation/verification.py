import copy, numpy as np, direct_integration as di, run_cases as rc
# 1) analytic: pseudo-first-order e- scavenging (only R12 active, Cr3+ in excess), instantaneous pulse
cfg = rc.case_config("Cr3+", 3.0); cfg["background_first_order_s^-1"]={}
cfg["disabled_reactions"]=["R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","R11"]
s=di.build_system(cfg); t=np.linspace(0,50e-9,2001); tt,Y,_=di.integrate_system(s,50e-9,t_eval=t)
C0=s.yields["e_s-"]; ce=Y[:,s.species_index["e_s-"]]; cr=Y[:,s.species_index["Cr3+"]]
# exact: d e/dt = -k e Cr, Cr = Cr0 - (C0 - e)  -> mixed second order A+B
k=6.1e10; a0=3e-3; b0=C0
exact = b0*(a0-b0)/(a0*np.exp((a0-b0)*k*tt)-b0)
print("analytic A+B (R12):   max rel err", np.max(np.abs(ce-exact))/C0)
# 2) analytic: pure second-order 2A->products (R5 only), Cl2•- only
cfg = rc.case_config("Cr3+", 0.0); cfg["background_first_order_s^-1"]={}
cfg["G_values_override"]={"Cl2•-":2.8}
cfg["disabled_reactions"]=[r for r in ["R1","R2","R3","R4","R6","R7","R8","R9","R10","R11","R12"]]
s=di.build_system(cfg); t=np.linspace(0,1e-3,2001); tt,Y,_=di.integrate_system(s,1e-3,t_eval=t)
C0=s.yields["Cl2•-"]; c=Y[:,s.species_index["Cl2•-"]]; exact=C0/(1+2*2.2e9*C0*tt)
print("analytic 2A (R5):     max rel err", np.max(np.abs(c-exact))/C0)
# 3) analytic: first-order background only
cfg = rc.case_config("Cr2+", 0.0); cfg["disabled_reactions"]=[f"R{i}" for i in range(1,13)]
s=di.build_system(cfg); t=np.linspace(0,100e-9,2001); tt,Y,_=di.integrate_system(s,100e-9,t_eval=t)
C0=s.yields["e_s-"]; c=Y[:,s.species_index["e_s-"]]; exact=C0*np.exp(-2.2e7*tt)
print("analytic 1st order:   max rel err", np.max(np.abs(c-exact))/C0)
# 4) tolerance convergence on full network, worst case 5 mM Cr3+ and 0.99 Cr2+
for metal,c0,T in [("Cr3+",5.0,300e-9),("Cr2+",0.99,50e-9)]:
    cfg=rc.case_config(metal,c0); s=di.build_system(cfg); t=np.linspace(0,T,3001)
    ref=di.integrate_system(s,T,t_eval=t,rtol=1e-12,atol=1e-18)[1]
    for rt,at in [(1e-6,1e-12),(1e-8,1e-14),(1e-10,1e-16)]:
        Y=di.integrate_system(s,T,t_eval=t,rtol=rt,atol=at)[1]
        e=[np.max(np.abs(Y[:,s.species_index[sp]]-ref[:,s.species_index[sp]]))/np.max(ref[:,s.species_index[sp]]) for sp in ("e_s-","Cl2•-")]
        print(f"conv {metal} {c0} rtol={rt:g}: max rel err e={e[0]:.1e} Cl2={e[1]:.1e}; min conc {Y.min():.1e}")
import pandas as pd
et,ea=rc.load_exp("absorbance5mMCr3.csv")
v=[rc.simulate_at(rc.case_config("Cr3+",5.0),et,n_out=n)[0] for n in (4000,20000,100000)]
print("grid: max rel diff 4k vs 100k", np.max(np.abs(v[0]-v[2]))/v[2].max(), " 20k vs 100k", np.max(np.abs(v[1]-v[2]))/v[2].max())
