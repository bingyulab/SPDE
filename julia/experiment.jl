## Model of Central Metabolic Pathway:  Glycolisis + TCA + ER Calcium crosstalk
# Author : Ingrid LIZCANO PRADA, Bingyu JIANG
# Last update : 10/27/2025
# import Pkg; Pkg.add("JumpProcesses")
using Random
using Timers
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DifferentialEquations
using StochasticDiffEq
using Statistics
using Plots
using DataFrames
using CSV
using StatsPlots
using Sundials

tic()

include("parameters.jl")

println("Building symbolic stochastic model...")

# 1. additive parameters with 0.5 strength ✅ unstable, too large
# 2. state-dependent and multiply parameters with 0.5 strength ✅
# 3. jump process ✅
# 4. modify ip3 from 0-2.
# 5. remove colored noise ✅
# 6. Caer, cac, atpc, adpc, psi, pyrm plot in the figures. ✅

# ============================================================================
# NOISE PARAMETERS (adjust these for different noise strengths)
# ============================================================================
noise_list = [:none, :additive, :multiplicative, :state_dependent, :jump]
colors = [:black, :blue, :red, :green, :purple, :orange]
important_variables = ["caer", "cac", "atpc", "adpc", "psi", "pyrm"]

# Gaussian Noises
const σ_additive =  0.05       # Additive Gaussian noise strength, 0.5, 0.1 are unstatble
const σ_multiplicative = 0.5   # Multiplicative Gaussian noise strength
const σ_calcium = 0.5          # State-dependent noise:Noise strength for calcium-dependent noise
const λ_jump = 0.01            # Jump rate (jumps per unit time)
const σ_jump = 0.01            # Jump size

# ============================================================================
# BUILD SYMBOLIC MODEL WITH NOISE OPTIONS
# ============================================================================

function build_stochastic_model(noise_type::Symbol)
    """
    Build symbolic stochastic differential equation model
    noise_type options:
    - :none          - Deterministic (no noise)
    - :additive      - Additive Gaussian noise on Jerout
    - :multiplicative - Multiplicative noise on Jerout  
    - :calcium       - State-dependent noise on calcium fluxes
    """
    
    # Define variables
    vars = @variables begin
        #Variable metabolite
        adpc(t)  
        adpm(t)  
        akg(t)  
        atpc(t)  
        atpm(t)
        cac(t)   
        cam(t)   
        fum(t)  
        isoc(t)  
        mal(t)
        nadhm(t) 
        nadm(t)  
        oaa(t)  
        psi(t)   
        scoa(t)
        suc(t)   
        cit(t)   
        caer(t)
        accoa(t)
        pyrm(t)
        pyrc(t)
        PEP(t)
        PG2(t)
        PG3(t)
        B13PG(t)
        GAP(t)
        DHAP(t)
        F16B(t)
        F6P(t)
        G6P(t)
        Gluc(t)

        # Algebraic helpers for flux calculation
        VDf1(t)
        Af1(t) 
        deltaMuH(t)
        VDres(t) 
        VDuni(t) 
        trc(t) 
        mwc(t) 
        Ares(t)
        Nhx(t)
        Dpfk2(t)
        Npfk(t)
        gapdh_denom(t)
        d_ldh(t)
        U(t)

        #Fluxes
        Jglctr(t)
        Jhx(t)
        Jgpi(t)
        Jpfk2(t)
        Jf26b(t)
        Jpfk(t)
        Jald(t)
        Jtpi(t)
        Jgapdh(t)
        Jpgk(t)
        Jpgm(t)
        Jeno(t)
        Jpk(t)
        Jldh(t)
        Jpyrh(t)
        Jpdhf(t)
        Jcs(t) 
        Jaco(t) 
        Jidh(t) 
        Jkgdh(t) 
        Jsl(t) 
        Jsdh(t) 
        Jfh(t) 
        Jmdh(t)
        Jf1(t) 
        Jhl(t) 
        Jhyd(t) 
        Jo(t) 
        Jant(t) 
        Jerout(t) 
        Jncx(t) 
        Jserca(t) 
        Juni(t)
        Jatpex(t)
        Jaat(t)

        #Flux with scales
        Jglctrs(t)
        Jhxs(t)
        Jgpis(t)
        Jpfk2s(t)
        Jf26bs(t)
        Jpfks(t)
        Jalds(t)
        Jtpis(t)
        Jgapdhs(t)
        Jpgks(t)
        Jpgms(t)
        Jenos(t)
        Jpks(t)
        Jldhs(t)
        Jpyrhs(t)
        Jf1s(t)
        Jants(t)
        Jhyds(t)
        Jpdhfs(t)
    end
    # Define ip3 as a symbolic parameter for optimization/scanning
    @parameters ip3
    @brownians B  # Definie Brownian Motion
    # ============================================================================
    # FLUX EQUATIONS (same as original)
    # ============================================================================
    
    flux_eqs = [
        #Glucose Transporter - Jglctr - (Berndt et al., 2015)
        Jglctr ~ ((V_mf_glut*(GLCex/k_glc_glut))-(V_mr_glut*(Gluc/k_glc_glut)))/(1.0 + (GLCex/k_glc_glut) + (Gluc/k_glc_glut))
        #Hexokinase - Jhx - (Mulukutla et al., 2015)
        Nhx ~ 1+(MgATP/k_i_mgatp_hk)+(G6P/k_i_g6p_hk)+(Gluc/k_m_glc_hk)+((MgATP*Gluc)/(k_m_mgatp_hk*k_m_glc_hk))+(MgADP/k_i_mgadp_hk) + (MgADP*G6P/(k_i_mgadp_hk*k_m_g6p_hk))+(Gluc*G6P/(k_m_glc_hk*k_m_g6p_hk)) + (Gluc*Ccg16p/(k_m_glc_hk*k_i_g16bp_hk)) + (Gluc*Cc23p2g/(k_m_glc_hk*k_i_23bpg_hk)) + (Gluc*Ccgsh/(k_m_glc_hk*k_i_gsh_hk))               
        Jhx ~ ((V_mf_hk*MgATP*Gluc/(k_m_mgatp_hk*k_m_glc_hk))-(V_mr_hk*MgADP*G6P/(k_i_mgadp_hk*k_m_g6p_hk)))*(1/Nhx)
        #Glucose-6-phosphate isomerase - Berndt 2015
        Jgpi ~ V_max_g6pi*((G6P-(F6P/k_e_g6pi))/(1+(G6P/k_m_g6p_g6pi)+(F6P/k_m_f6p_g6pi)))
        # 6-Phosphofructo-2-Kinase - Jpfk2 -
        Dpfk2 ~ (k_i_atp_pfk2*k_m_f6p_pfk2)+(k_m_f6p_pfk2*atpc)+(k_m_atp_pfk2*F6P)+(k_m_adp_pfk2*F26B/k_eq_pfk2)+ (k_m_f26p_pfk2*adpc/k_eq_pfk2)+ (atpc*F6P)+ (k_m_adp_pfk2*atpc*F26B/(k_eq_pfk2*k_i_atp_pfk2))+ (adpc*F26B/(k_eq_pfk2)) + (k_m_atp_pfk2*adpc*F26B/k_i_adp_pfk2)+ (atpc*F6P*F26B/k_i_f26p_pfk2) + (adpc*F6P*F26B/(k_eq_pfk2*k_i_f6p_pfk2))
        Jpfk2 ~ (V_f_pfk2*((atpc*F6P)-(adpc*F26B/k_eq_pfk2)))/(Dpfk2*(1+(PEP/k_i_pep_pfk2)))
        #Fructose-2,6-Bisphosphatase - Jf26b -
        Jf26b ~ (V_m_f26bp*F26B)/((1+(F6P/k_i_f6p_f26b))*(k_m_f26p_f26b+ F26B))
        #Phosphofructokinase - Jpfk - (Mulukutla et al., 2015)
        Npfk ~ 1 + ((Lpfk*((1+(atpc/k_m_atp_pfk))^4))*((1+(Mg/k_m_mg_pfk))^4)*((1+(Cc23p2g/k_m_23bpg_pfk))^4)) / (((1+(F6P/k_m_f6p_pfk)+(F16B/k_m_f16p_pfk))^4)*((1+(Ccamp/k_m_amp_pfk))^4)*((1+(Ccg16p/k_m_g16bp_pfk))^4)*((1+(Ccpi/k_m_pi_pfk))^4)*((1+(F26B/k_m_f26p_pfk))^4))  
        Jpfk ~ (((V_f_pfk*MgATP*F6P/(k_m_f6p_pfk*k_m_mgatp_pfk))-(V_r_pfk*MgADP*F16B/(k_m_f16p_pfk*k_m_mgadp_pfk)))/(((1+(F6P/k_m_f6p_pfk))*(1+(MgATP/k_m_mgatp_pfk)))+((1+(F16B/k_m_f16p_pfk))*(1+(MgADP/k_m_mgadp_pfk)))-1)) *(1/Npfk)
        #Aldolase - Jald -
        Jald ~ ((V_mf_ald*F16B/k_m_f16p_ald)-(V_mr_ald*GAP*DHAP/(k_m_gap_ald*k_i_dhap_ald)))/ (1+(Cc23p2g/k_i_23bpg_ald)+(F16B/k_m_f16p_ald)+((k_m_dhap_ald*GAP/(k_m_gap_ald*k_i_dhap_ald))*(1+(Cc23p2g/k_i_23bpg_ald)))+(DHAP/k_i_dhap_ald)+(k_m_dhap_ald*F16B*GAP/(k_i_f16p_ald*k_m_gap_ald*k_i_dhap_ald))+(DHAP*GAP/(k_m_gap_ald*k_i_dhap_ald)))
        #Triose Phosphate Isomerase - Jtpi -
        Jtpi ~ ((V_mf_tpi*DHAP/k_f_tpi)-(V_mr_tpi*GAP/k_r_tpi))/(1+(DHAP/k_f_tpi)+(GAP/k_r_tpi))
        #Glyceraldehyde 3-phosphate dehydrogenase - Jgapdh - berndt
        gapdh_denom ~ ((1.0+(nadc/k_m_nad_gapdh))*(1.0+(GAP/k_m_grap_gapdh))*(1.0+(Ccpi/k_m_pic_gapdh)))+((1.0+(nadhc/k_m_nadh_gapdh))*(1.0+(B13PG/k_m_bpg13_gapdh)))-1.0
        Jgapdh ~ V_max_gapdh*(((nadc*GAP*Ccpi)-(B13PG*nadhc/k_e_gapdh))/(gapdh_denom))
        #Phosphoglycerate Kinase - Jpgk -
        Jpgk ~ ((V_mf_pgk*B13PG*MgADP/(k_i_mgadp_pgk*k_m_b13_pgk))-(V_mr_pgk*PG3*MgATP/(k_i_mgatp_pgk*k_m_pg3_pgk)))/(1+(B13PG/k_i_b13_pgk)+(MgADP/k_i_mgadp_pgk)+(B13PG*MgADP/(k_i_mgadp_pgk*k_m_b13_pgk))+(PG3/k_i_pg3_pgk)+(MgATP/k_i_mgatp_pgk)+(PG3*MgATP/(k_i_mgatp_pgk*k_m_pg3_pgk)))
        #Phosphoglycerate Mutase - Jpgm -
        Jpgm ~ ((V_mf_pgm*PG3/k_m_pg3_pgm)-(V_mr_pgm*PG2/k_m_pg2_pgm))/(1+(PG3/k_m_pg3_pgm)+(PG2/k_m_pg2_pgm))
        #Enolase - Jeno -
        Jeno ~ ((V_mf_eno*PG2*Mg/(k_i_mg_eno*k_m_pg2_eno))-(V_mr_eno*PEP*Mg/(k_i_mg_eno*k_m_pep_eno)))/ (1+ (PG2/k_i_pg2_eno)+(Mg/k_i_mg_eno)+(PG2*Mg/(k_i_mg_eno*k_m_pg2_eno))+(PEP/k_i_pep_eno)+(Mg/k_i_mg_eno)+(PEP*Mg/(k_i_mg_eno*k_m_pep_eno)))
        #Pyruvate kinase - Jpk - Berndt 2015
        Jpk ~ V_max_pkb*(PEP/(PEP+k_m_pep_pkb))*(adpc/(adpc+k_m_adp_pkb*(1+(atpc/k_i_atp_pkb))))
        # Lactate Dehydrogenase- Jldh - 
        d_ldh ~ (nadhc/k_i_nadh_ldh)+(nadc/k_i_nad_ldh)+((nadhc*pyrc)/(k_i_nadh_ldh*k_m_pyr_ldh))+((k_m_nad_ldh*nadhc*LAC)/(k_i_nad_ldh*k_i_nadh_ldh*k_m_lac_ldh))+ ((k_m_nadh_ldh*nadc*pyrc)/(k_i_nad_ldh*k_i_nadh_ldh*k_m_pyr_ldh))+((nadc*LAC)/(k_i_nad_ldh*k_m_lac_ldh))+((nadhc*pyrc*LAC)/(k_i_nadh_ldh*k_m_pyr_ldh*k_i_lac_ldh))+((nadc*pyrc*LAC)/(k_i_nad_ldh*k_i_pyr_ldh*k_m_lac_ldh))
        Jldh ~ ((V_mf_ldh*nadhc*pyrc/(k_i_nadh_ldh*k_m_pyr_ldh))-(V_mr_ldh*nadc*LAC/(k_i_nad_ldh*k_m_lac_ldh)))/(((1+((k_m_nadh_ldh*pyrc)/(k_i_nadh_ldh*k_m_pyr_ldh))+((k_m_nad_ldh*LAC)/(k_i_nad_ldh*k_m_lac_ldh)))*(1+(pyrc/k_i_pyr_ldh)))+d_ldh)
        #Pyruvate –Hydrogen shuttle - Jpyrh - 
        Jpyrh ~ V_max_pyrh*((pyrc*hc)-(pyrm*hm))
        # Pyruvatedehydrogenase complex - Jpdh - 
        Jpdhf ~ V_max_pdhf*(1+(A_max_cam_pdhf*cam/(cam+k_a_cam_pdhf)))*(pyrm/(pyrm+k_m_pyr_pdhf))*(fadm/(fadm+k_m_NADm_pdhf))*(coa/(coa+(k_m_coam_pdhf*(1.0 +(accoa/k_i_accoa_pdhf)))))     # (Berndt et al., 2015)
        
        # TCA / fluxes
        Jcs  ~ Vmaxcs*accoa/(accoa+kmaccoa+((kmoaa*accoa/oaa)*(1+accoa/kiaccoa))+(ksaccoa*kmoaa)/oaa)
        Jaco ~ kfaco*(cit-isoc/keaco)
        Jidh ~ Vmaxidh/(1+(hm/kh1)+(kh2/hm)+((kmisoc/isoc)^ni)/((1+adpm/kaadp)*(1+cam/kaca))
                        +(kmnad/nadm)*(1+nadhm/kinadh)
                        +(((kmisoc/isoc)^ni)*(kmnad/nadm)*(1+nadhm/kinadh))/((1+adpm/kaadp)*(1+cam/kaca)))
        
        Jkgdh ~ Vmaxkgdh/(1+((kmakg/akg)*((kmnadkgdh/nadm)^nakg))/((1+Mg/kdmg)*(1+cam/kdca)))
        Jsl ~ kfsl*(scoa*adpm*pim-suc*atpm*coa/kesl)
        Jsdh ~ Vmaxsdh/(1+(kmsuc/suc)*(1+oaa/kioaasdh)*(1+fum/kifum))
        Jfh ~ kffh*(fum-mal/kefh)
        Jmdh ~ Vmdh*(mal*nadm-(oaa*nadhm/Keqmdh))/((1+mal/kmmal)*(1+nadm/kmnadmdh)+(1+oaa/kmoaamdh)*(1+nadhm/kmnadhmdh)-1)
        
        # Algebraic constants for Phosphorylation of mitochondrial ADP via F1F0 ATPase
        VDf1  ~ exp((3.0*psi)/RTF)
        Af1  ~ kf1*atpm/(adpm*(pim/uMmM))         # pim in mM; uMmM converts uM terms
        deltaMuH ~ psi-(2.303*RTF*DpH)
        VDres    ~ exp(6.0*g*psi/RTF)
        
        # F1F0-ATPase, leak, hydrolysis, respiration, ANT
        Jf1  ~ -rhof1*((Pac1*Af1)-(pa*VDf1)+(pc2*Af1*VDf1))/(((1.0+p1*Af1)*VDf1B)+((p2+p3*Af1)*VDf1))
        Jhl  ~ gH*deltaMuH
        Jhyd ~ khyd*atpc/(katpcH + atpc)
        
        # Oxidation of NADH by respiration
        Ares ~ kres*(nadhm^0.5/nadm^0.5) 
        Jo   ~ 0.5 * rhores * ((Rac1 * Ares) - (ra * VDres) + (rc2 * Ares * VDres)) /
            ( ((1.0 + r1 * Ares) * VDresB) + ((r2 + r3 * Ares) * VDres) )
        Jant ~ vant * (1.0 - alphac * atpc * adpm * exp(-psi / RTF) / (alpham * adpc * atpm)) /
            ( (1.0 + alphac * atpc * exp(-fant * psi / RTF) / adpc) * (1. + adpm / (alpham * atpm)) )
        
        # ER/IP3R (TO ADD NOISE IN Jerout)
        Jerout ~ ((vip3 * ((ip3^2.0)/((ip3^2.0)+(kip3^2.0))) * ((cac^2.0)/((cac^2.0)+(kcaa^2.0))) *
                ((kcai1^4.0)/((kcai2^4.0)+(cac^4.0)))) + vleak) * (caer - cac) / uMmM
        
        # NCX + SERCA 
        Jncx   ~ vncx * (cam / cac) * exp(b * (psi - psiStar) / RTF) / ((1. + kna / nac)^n * (1. + kca / cam))
        Jserca ~ vserca * (cac^2. / (kserca^2. + cac^2.)) * (atpc / (atpc + katpc))
        
        # MCU helpers
        VDuni ~ 2.0 * (psi - psiStar) / RTF
        trc ~ cac / ktrans
        mwc ~ trc * ((1+ trc)^3.0) / (((1.0 + trc)^4.0) + (L / ((1.0 + (cac / kact))^ma)))
        Juni ~ (vuni * VDuni * mwc) / (1.0 - exp(-VDuni))

        # ATP-exchanger    #in is assumed as cytosol
        U ~ psi/(1000*RTF)
        Jatpex ~ V_max_atpex*((1.0-(atpc*adpm/(adpc*atpm))*exp(U))/(1.0+((atpc/adpc)*exp(Spsi*U)*(1.0+(adpm/atpm)))))

        # Mitochondrial aspartate aminotransferase
        Jaat ~ V_max_aatm*((aspm*akg-(oaa*glutm/k_e_aatm)))
        
        #scales
        Jf1s ~ s_f1 * Jf1
        Jants ~ s_ant * Jant
        Jhyds ~ s_hyd * Jhyd
        Jpdhfs ~ s_pdh * Jpdhf
        Jglctrs ~ s_glctr * Jglctr
        Jhxs ~ s_hx * Jhx
        Jgpis ~ s_gpi * Jgpi
        Jpfk2s ~ s_pfk2 * Jpfk2
        Jf26bs ~ s_f26b * Jf26b
        Jpfks ~ s_pfk * Jpfk
        Jalds ~ s_ald * Jald
        Jtpis ~ s_tpi * Jtpi
        Jgapdhs ~ s_gapdh * Jgapdh
        Jpgks ~ s_pgk * Jpgk
        Jpgms ~ s_pgm * Jpgm
        Jenos ~ s_eno * Jeno
        Jpks ~ s_pk  * Jpk
        Jldhs ~ s_ldh  * Jldh
        Jpyrhs ~ s_pyrh  * Jpyrh
    ]
    # ============================================================================
    # BUILD NOISE EQUATIONS BASED ON TYPE
    # ============================================================================
          
    d_caer_base = uMmM*deltaer*(Jserca - Jerout)
    d_cac_base = uMmM*fc*(-Jserca + Jerout + delta * (Jncx - Juni))  

    if noise_type in [:none, :jump] 
        d_caer_eq = D(caer) ~ d_caer_base
        d_cac_eq = D(cac) ~ d_cac_base
    # elseif noise_type == :colored
    #     @variables η(t)        
    #     d_caer_eq = D(caer) ~ d_caer_base + (- k_ou * η)
    #     d_cac_eq = D(cac) ~ d_cac_base + (k_ou * η)
    elseif noise_type == :additive
        d_caer_eq = D(caer) ~ d_caer_base + (-σ_additive * B)
        d_cac_eq = D(cac) ~ d_cac_base + (σ_additive * B)

    elseif noise_type == :multiplicative

        d_caer_eq = D(caer) ~ d_caer_base + (-σ_multiplicative * Jerout * B)
        d_cac_eq = D(cac) ~ d_cac_base + (σ_multiplicative * Jerout * B)
    elseif noise_type == :state_dependent
        d_caer_eq = D(caer) ~ d_caer_base + (-σ_calcium * sqrt(abs(Jerout)) * B)
        d_cac_eq = D(cac) ~ d_cac_base + (σ_calcium * sqrt(abs(Jerout)) * B)
    else        
        error("Unsupported noise type: $noise_type")
    end

    diff_eqs = [
            # Differential equations for metabolite concentration
            D(adpc) ~ -delta * Jants + Jhyds + 0.5*Jserca                       #ADP_c mM/s
            D(adpm) ~ Jants - Jf1s - Jsl                                        # adpm mM/s
            D(akg) ~ Jidh - Jkgdh #+ Jaat                                     # akg mM/s
            D(atpc) ~ -(-delta * Jants + Jhyds + 0.5*Jserca)                    #ATP_c mM/s
            D(atpm) ~ -(Jants - Jf1s - Jsl)                                     # atpm mM/s
            d_cac_eq                                                           # cac uM/s
            D(cam) ~ uMmM*fm*(Juni - Jncx)                                    # cam uM/s
            D(cit) ~ Jcs - Jaco                                               # cit mM/s
            D(fum) ~ Jsdh - Jfh                                               # fum mM/s
            D(isoc) ~ Jaco - Jidh                                             # isoc mM/s
            D(mal) ~  Jfh - Jmdh                                              # mal mM/s
            D(nadm) ~ -(-Jo + Jidh + Jkgdh + Jmdh)                            # nadm mM/s
            D(nadhm) ~ -Jo + Jidh + Jkgdh + Jmdh                              # nadhm mM/s
            D(oaa) ~ Jmdh - Jcs #- Jaat                                       # oaa mM/s
            D(psi) ~ (10*Jo - 3*Jf1s - Jants - Jhl - Jncx - 2*Juni)/cmito       # psi mV/s
            D(scoa) ~ Jkgdh - Jsl                                             # scoa mM/s
            D(suc) ~ Jsl - Jsdh                                               # suc mM/s
            d_caer_eq                                                           # Caer uM/s
            D(accoa) ~ Jpdhfs - Jcs                                            #accoa mM/s
            D(pyrm) ~ Jpyrh - Jpdhfs                                           #pyrm mM/s
            D(pyrc) ~ Jpks - Jldh - Jpyrh                                      #Pyrc
            D(PEP) ~ Jeno - Jpks                                               #PEP
            D(PG2) ~ Jpgm - Jeno                                              #PG2
            D(PG3) ~ Jpgks - Jpgm                                              #PG3
            D(B13PG) ~ Jgapdh - Jpgks                                          #B13PG
            D(GAP) ~ Jald + Jtpi - Jgapdh                                     #GAP
            D(DHAP) ~ Jald - Jtpi                                             #DHAP
            D(F16B) ~ Jpfks - Jald                                             #F16B
            D(F6P) ~ Jgpi  - Jpfks  - Jpfk2 + Jf26b                            #F6P
            D(G6P) ~ Jhxs - Jgpi                                               #G6P
            D(Gluc) ~ Jglctrs - Jhxs                                            #Gluc
        ]    
        
    # if noise_type in [:colored]        
    #     push!(diff_eqs, D(η) ~ -η/τ_ou + σ_ou * B)  
    # end

    # ============================================================================
    # INITIAL CONDITIONS (moved here to be in scope of symbolic variables)
    # ============================================================================
    
    adpc0 = atot*0.20
    adpm0 = amtot*0.50 
    akg0  = ckint*0.01 
    atpc0 = atot - adpc0
    atpm0 = amtot - adpm
    cac0  = 0.2
    cam0  = 0.10
    fum0  = ckint*0.01
    isoc0 = ckint*0.01
    mal0  = ckint*0.01
    nadhm0 = 0.125*nadtot
    nadm0  = nadtot - 0.125*nadtot
    oaa0  = 0.001
    psi0  = 160.0
    scoa0 = ckint*0.01
    suc0  = ckint*0.01
    cit0  = ckint - isoc0 - akg0 - scoa0 - suc0 - fum0 - mal0 - oaa0
    caer0 = coef1 - cac0 * coef2 - cam0 * coef3
    accoa0= 0.01
    pyrm0= 5 
    pyrc0 =5
    PEP0 = 0.017 #Clement 2020 #2.67 #Chassgnole et al. 2002 (mM) measured
    PG20 = 0.01 #Clement 2020 #0.399 #Chassgnole et al. 2002 (mM) measured
    PG30 = 0.069 #2.13 #Chassgnole et al. 2002 (mM) measured
    B13PG0 = 0.000369 #0.008 #Chassgnole et al. 2002 (mM) measured
    GAP0 = 0.00194 #0.048 #Penkler et al., 2014
    DHAP0 = 0.02
    F16B0 = 0.00231 #1.0 #0.272
    F6P0 = 0.013 #Clement 2020 #2.0 #Penkler et al., 2014 0.24 #Penkler et al., 2014
    G6P0 = 0.039  #Clement 2020 #G6P0=3.48
    Gluc0 = 10.0 #Clement 2020
    
    u0 = Dict(
        adpc => adpc0,  adpm => adpm0,  akg => akg0,   atpc => atpc0, atpm => atpm0,
        cac => cac0,    cam => cam0,    fum => fum0,   isoc => isoc0, mal => mal0,
        nadhm => nadhm0, nadm => nadm0, oaa => oaa0,   psi => psi0,   scoa => scoa0,
        suc => suc0,    cit => cit0,    caer => caer0,    accoa => accoa0,    pyrm => pyrm0,
        pyrc => pyrc0,    PEP => PEP0,    PG2 => PG20,    PG3 => PG30,    B13PG => B13PG0,
        GAP => GAP0,    DHAP => DHAP0,    F16B => F16B0,    F6P => F6P0,    G6P => G6P0,   Gluc => Gluc0
    )

    # if noise_type == :colored
    #     u0[η] = 0.0
    # end

    # Combine all equations
    all_eqs = vcat(flux_eqs, diff_eqs)

    return all_eqs, u0, ip3
end

const var_names = ["adpc", "adpm", "akg", "atpc", "atpm",
                   "cac", "cam", "fum", "isoc", "mal",
                   "nadhm", "nadm", "oaa", "psi", "scoa",
                   "suc", "cit", "caer", "accoa", "pyrm",
                   "pyrc", "PEP", "PG2", "PG3", "B13PG",
                   "GAP", "DHAP", "F16B", "F6P", "G6P", "Gluc"]


function simulate_model(noise_type::Symbol; tspan=(0.0, 4000.0), ip3_val=0.7)
    println("\n" * "="^80)
    println("SIMULATING MODEL WITH NOISE TYPE: $noise_type")
    println("="^80)
    
    # CRITICAL FIX: Use different random seed for each noise type
    seeds = Dict(
        :none => 1111, :additive => 1234, :multiplicative => 5678,
        :state_dependent => 9012, :colored => 3456
    )
    current_seed = get(seeds, noise_type, 1111)
    Random.seed!(current_seed)
    println("USING RANDOM SEED: $current_seed")

    # Build symbolic model
    all_eqs, u0, ip3_param = build_stochastic_model(noise_type)

    # Create system
    if noise_type in [:none, :jump] # ,
        println("Solving ODE system...") 
        @named sys_raw = ODESystem(all_eqs, t)
        sys_raw = structural_simplify(sys_raw)  
    elseif noise_type in [:additive, :multiplicative, :state_dependent]
        println("Using brownian motion...")
        @mtkcompile sys_raw = System(all_eqs, t)
    else
        error("Unknown noise type: $noise_type")
    end
    
    caer_idx = findfirst(eq -> occursin("Differential(t)(caer", string(eq.lhs)), equations(sys_raw))
    cac_idx = findfirst(eq -> occursin("Differential(t)(cac", string(eq.lhs)), equations(sys_raw))            
    println("DEBUG $noise_type: Caer with index $caer_idx is defined as $(equations(sys_raw)[caer_idx])")
    println("DEBUG $noise_type: Cac with index $cac_idx is defined as $(equations(sys_raw)[cac_idx])")
    
    u0 = merge(u0, Dict(ip3_param => ip3_val))

    # Solve
    if noise_type == :none
        println("Solving deterministic ODE...")
        prob = ODEProblem(sys_raw, u0, tspan)
        sol = solve(prob; reltol=1e-6, abstol=1e-9)
    elseif noise_type == :additive 
        prob = SDEProblem(sys_raw, u0, tspan)
        sol = solve(prob, ImplicitRKMil(), adaptive=true, saveat=1.0,
                    reltol=1e-4, abstol=1e-6, maxiters=10^7,
                    seed=current_seed,
                    isoutofdomain=(u,p,t) -> any(isnan.(u)) || any(u .< 0))
        
    elseif noise_type == :multiplicative 
        prob = SDEProblem(sys_raw, u0, tspan)
        sol = solve(prob, ImplicitRKMil(), adaptive=true, saveat=1.0,
                    reltol=1e-4, abstol=1e-6, maxiters=10^7,
                    seed=current_seed,
                    isoutofdomain=(u,p,t) -> any(isnan.(u)) || any(u .< 0))
    
    elseif noise_type == :state_dependent 
        
        prob = SDEProblem(sys_raw, u0, tspan)
        sol = solve(prob, ImplicitRKMil(), adaptive=true, saveat=1.0,
                    reltol=1e-4, abstol=1e-6, maxiters=10^7,
                    seed=current_seed,
                    isoutofdomain=(u,p,t) -> any(isnan.(u)) || any(u .< 0))
    elseif noise_type == :colored
        prob = SDEProblem(sys_raw, u0, tspan)
        sol = solve(prob, ISSEM(), adaptive=true, saveat=1.0,
                    reltol=1e-5, abstol=1e-8, maxiters=10^7,
                    seed=current_seed,
                    isoutofdomain=(u,p,t) -> any(isnan.(u)) || any(u .< 0))
    elseif noise_type == :jump
        println("Solving Jump SDE...")

        prob = ODEProblem(sys_raw, u0, tspan)
        # prob = SDEProblem(sys_raw, u0, tspan)
        # Define jump process
        rate1(u, p, t) = u[caer_idx] * λ_jump  
        rate2(u, p, t) = u[cac_idx] * λ_jump  
        function affect!(integrator) 
            # Add random jump to caer
            integrator.u[caer_idx]  += σ_jump
            integrator.u[cac_idx] -= σ_jump
            nothing
        end

        jump1 = ConstantRateJump(rate1, affect!)
        jump2 = ConstantRateJump(rate2, affect!)
        jump_prob = JumpProblem(prob, Direct(), jump1, jump2)
        sol = solve(jump_prob, Rodas5P(); saveat=1.0, 
                    isoutofdomain=(u,p,t) -> any(isnan.(u)) || any(u .< 0))
    else        
        error("Unknown noise type: $noise_type")
    end
    
    println("Simulation complete!")
    return sol, var_names
end

# ============================================================================
# ADDITIONAL DIAGNOSTIC: Verify noise is actually different
# ============================================================================

function verify_noise_differences(results)
    """
    More stringent verification that noise types produce different results
    """
    println("\n" * "="^80)
    println("VERIFICATION: Are noise types producing different trajectories?")
    println("="^80)
    
    # Extract last 500 points of caer for each noise type
    trajectories = Dict()
    
    for (nt, (sol, _)) in results
        df = DataFrame(sol)
        if nrow(df) < 500
            continue
        end
        
        caer_idx = df_find_column(df, "caer")
        if isnothing(caer_idx)
            continue
        end
        
        # Take last 500 points
        caer_data = df[end-499:end, caer_idx]
        trajectories[nt] = caer_data
    end
    
    # Compare each pair
    pairs_tested = 0
    pairs_different = 0
    
    noise_types = collect(keys(trajectories))
    
    for i in 1:length(noise_types)
        for j in (i+1):length(noise_types)
            nt1, nt2 = noise_types[i], noise_types[j]
            
            data1 = trajectories[nt1]
            data2 = trajectories[nt2]
            
            # Statistical tests
            correlation = cor(data1, data2)
            mae = mean(abs.(data1 .- data2))
            max_diff = maximum(abs.(data1 .- data2))
            
            # Kolmogorov-Smirnov test would be ideal here
            # For now, use simple metrics
            
            pairs_tested += 1
            
            println("\n$nt1 vs $nt2:")
            println("  Correlation: $(round(correlation, digits=4))")
            println("  MAE: $(round(mae, digits=6))")
            println("  Max difference: $(round(max_diff, digits=6))")
            
            # Consider different if correlation < 0.95 OR MAE > 0.1
            if correlation < 0.95 || mae > 0.1
                pairs_different += 1
                println("  ✅ DIFFERENT")
            else
                println("  ⚠️  SIMILAR")
            end
        end
    end
    
    println("\n" * "-"^80)
    println("Summary: $pairs_different/$pairs_tested pairs are detectably different")
    
    if pairs_different < pairs_tested / 2
        println("⚠️  WARNING: Many noise types appear too similar!")
        println("   Consider increasing noise strengths or checking implementation")
    else
        println("✅ Noise types are producing distinct trajectories")
    end
    
    println("="^80)
end


# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================
# Helper: find dataframe column index by matching name substring (case-insensitive)
function df_find_column(df, varname)
    for (i, nm) in enumerate(names(df))
        if occursin(varname, lowercase(string(nm)))
            return i
        end
    end
    return nothing
end

function plot_results(sol, var_names, noise_type)
    println("Creating plots...")
    dsol= DataFrame(sol)
    dsol_window = filter(:timestamp => t -> t ≥ (maximum(dsol.timestamp) - 1000), dsol)

    nvars = ncol(dsol_window) - 1  # exclude time column
    nrows = ceil(Int, sqrt(nvars))
    ncols = ceil(Int, nvars / nrows)
    
    plots = [
        plot(dsol_window.timestamp, dsol_window[!, i+1],
            xlabel="time", ylabel=names(dsol_window)[i+1], legend=false)
        for i in 1:nvars 
        if !(names(dsol_window)[i+1] in ["η(t)", "dummy_A(t)", "dummy_M1(t)", "dummy_M2(t)", "dummy_S1(t)", "dummy_S2(t)", "dummy_S3(t)"])
    ]

    p_detailed = plot(plots..., layout=(nrows, ncols), size=(1500,800))
    
    # Save
    mkpath("imgs/bio/symbolic/")
    savefig(p_detailed, "imgs/bio/symbolic/all_variables_$(noise_type).png")
    println("Saved: imgs/bio/symbolic/all_variables_$(noise_type).png")
    
    # Plot key variables
    key_indices = [df_find_column(dsol_window, v) for v in important_variables]

    plots_key = []
    for idx in key_indices
        if !isnothing(idx)
            col_idx = idx
            vname = string(names(dsol_window)[col_idx])
            p = plot(dsol_window.timestamp, dsol_window[!, col_idx],
                    xlabel="Time (s)", ylabel=vname, legend=false,
                    linewidth=2, title=vname)
            push!(plots_key, p)
        end
    end
    
    p_key = plot(plots_key..., layout=(2, 3), size=(1200, 600))
    savefig(p_key, "imgs/bio/symbolic/key_variables_$(noise_type).png")
    println("Saved: imgs/bio/symbolic/key_variables_$(noise_type).png")
    
    # Cac vs ATPc dual axis
    tg = dsol_window.timestamp    
    cac_idx = df_find_column(dsol_window, "cac")
    atpc_idx = df_find_column(dsol_window, "atpc")

    if !isnothing(cac_idx) && !isnothing(atpc_idx)
        cacg = dsol_window[!, cac_idx]
        atpg = dsol_window[!, atpc_idx]

        p_dual = plot(tg, cacg, color=:red, ylabel="Cac (mM)", 
                     xlabel="Time (s)", legend=false, linewidth=2)
        plot!(twinx(), tg, atpg, color=:blue, ylabel="ATPc (mM)", 
              legend=false, linewidth=2, title="Cac vs ATPc - $noise_type")
        
        savefig(p_dual, "imgs/bio/symbolic/cac_vs_atpc_$(noise_type).png")
        println("Saved: imgs/bio/symbolic/cac_vs_atpc_$(noise_type).png")
    end
    
    # ATPc/ADPc ratio
    tg = dsol_window.timestamp
    atpc_col = df_find_column(dsol_window, "atpc")
    adpc_col = df_find_column(dsol_window, "adpc")
    
    if !isnothing(atpc_col) && !isnothing(adpc_col)
        atpg = dsol_window[!, atpc_col]
        adpg = dsol_window[!, adpc_col]
        ratio = atpg ./ adpg
        
        # Filter out invalid values
        valid_idx = (atpg .> 0) .& (adpg .> 0) .& isfinite.(atpg) .& isfinite.(adpg)
        if sum(valid_idx) > 0
            ratio = atpg[valid_idx] ./ adpg[valid_idx]
            tg_valid = tg[valid_idx]
            
            p_ratio = plot(tg_valid, ratio, color=:purple, xlabel="Time (s)", 
                          ylabel="ATPc/ADPc", legend=false, linewidth=2,
                          title="ATP:ADP Ratio - $noise_type")
            
            savefig(p_ratio, "imgs/bio/symbolic/atpc_adpc_ratio_$(noise_type).png")
            println("Saved: imgs/bio/symbolic/atpc_adpc_ratio_$(noise_type).png")
        else
            println("Warning: No valid data for ATP:ADP ratio plot")
        end
    else
        println("Warning: Could not find ATP or ADP columns for ratio plot")
    end
    
    return p_detailed, p_key
end

function save_results(sol, var_names, noise_type)
    println("Saving results to CSV...")
    mkpath("results/symbolic/")
    
    df = DataFrame(sol)
    
    # Rename columns properly
    rename_dict = Dict(:timestamp => :time)
    for (i, vname) in enumerate(var_names)
        old_name = Symbol("value$i")
        if old_name in names(df)
            rename_dict[old_name] = Symbol(vname)
        end
    end
    rename!(df, rename_dict)
    
    CSV.write("results/symbolic/$(noise_type).csv", df)
    println("Saved: results/symbolic/$(noise_type).csv")
end

function compare_noise_types(results, var_names)
    """
    Create comparison plots for all noise types
    """
    println("\nCreating comparison plots...")
    
    # Focus on caer and cac
    caer_name = "caer"
    cac_name  = "cac"   
    
    # Create subplots for each noise type    
    plots_caer = []
    plots_cac = []
    
    for nt in noise_list
        if nt == :none
            continue  # Skip deterministic for this comparison
        end
        if haskey(results, nt)
            sol, _ = results[nt]
            dsol = DataFrame(sol)
            
            # Take last 1000s for steady-state comparison
            dsol_window = filter(:timestamp => t -> t ≥ (maximum(dsol.timestamp) - 1000), dsol)
            
            t = dsol_window.timestamp
            # SAFE: find columns by name in the dataframe rather than assuming fixed indices
            ccol = df_find_column(dsol_window, caer_name)
            kcol = df_find_column(dsol_window, cac_name)
            if isnothing(ccol) || isnothing(kcol)
                println("Warning: could not find caer/cac columns for $nt, skipping plot")
                continue
            end
            caer = dsol_window[!, ccol]
            cac  = dsol_window[!, kcol]
            
            p_caer = plot(t, caer, label=string(nt), 
                         xlabel="Time (s)", ylabel="CaER (μM)",
                         title="CaER - $nt", linewidth=1.5,
                         legend=false)
            
            p_cac = plot(t, cac, label=string(nt),
                        xlabel="Time (s)", ylabel="CaC (μM)", 
                        title="CaC - $nt", linewidth=1.5,
                        legend=false)
            
            push!(plots_caer, p_caer)
            push!(plots_cac, p_cac)
        end
    end
    cnt = ceil(Int, (length(noise_list)-1)/2)
    p_caer_combined = plot(plots_caer..., layout=(cnt, 2), size=(1400, 900))
    p_cac_combined = plot(plots_cac..., layout=(cnt, 2), size=(1400, 900))
    
    mkpath("imgs/bio/symbolic/")
    savefig(p_caer_combined, "imgs/bio/symbolic/comparison_caer_all_noise.png")
    savefig(p_cac_combined, "imgs/bio/symbolic/comparison_cac_all_noise.png")
    
    println("Saved comparison plots!")
    
    # Also create overlaid plot for direct comparison
    p_overlay_caer = plot(xlabel="Time (s)", ylabel="CaER (μM)", 
                          title="CaER Comparison - All Noise Types",
                          legend=:topright, size=(1000, 600))
    p_overlay_cac = plot(xlabel="Time (s)", ylabel="CaC (μM)",
                         title="CaC Comparison - All Noise Types", 
                         legend=:topright, size=(1000, 600))
        
    for (i, nt) in enumerate(noise_list)
        if haskey(results, nt)
            sol, _ = results[nt]
            dsol = DataFrame(sol)
            dsol_window = filter(:timestamp => t -> t ≥ (maximum(dsol.timestamp) - 1000), dsol)
            
            t = dsol_window.timestamp
            ccol = df_find_column(dsol_window, caer_name)
            kcol = df_find_column(dsol_window, cac_name)
            if isnothing(ccol) || isnothing(kcol)
                println("Warning: could not find caer/cac columns for $nt, skipping plot")
                continue
            end
            caer = dsol_window[!, ccol]
            cac  = dsol_window[!, kcol]
            
            plot!(p_overlay_caer, t, caer, label=string(nt), 
                  color=colors[i], linewidth=2, alpha=0.7)
            plot!(p_overlay_cac, t, cac, label=string(nt),
                  color=colors[i], linewidth=2, alpha=0.7)
        end
    end
    
    savefig(p_overlay_caer, "imgs/bio/symbolic/overlay_caer_comparison.png")
    savefig(p_overlay_cac, "imgs/bio/symbolic/overlay_cac_comparison.png")
    
    println("Saved overlay comparison plots!")
end


# ============================================================================
# MAIN EXECUTION
# ============================================================================

println("\n" * "="^80)
println("SYMBOLIC STOCHASTIC METABOLIC MODEL")
println("="^80 * "\n")

println("\n" * "="^80)
println("Using Brownian Motion Noise Implementation")
println("="^80)

results = Dict()

for noise_type in noise_list
    # Time the solve
    t_start = time()
    
    # Solve with current IP3 value
    sol, vnames = simulate_model(noise_type; tspan=(0.0, 4000.0))    
    t_solve = time() - t_start
    
    # Check convergence
    converged = (Symbol(sol.retcode) == :Success || Symbol(sol.retcode) == :Default) 
    final_t = sol.t[end]
    n_steps = length(sol.t)
    
    status_icon = converged ? "✓" : "⚠️"
    println("$status_icon Simulation $(sol.retcode): $(round(t_solve, digits=2))s, $n_steps steps, final t=$(round(final_t, digits=1))")
    
    results[noise_type] = (sol, vnames)
    
    plot_results(sol, vnames, noise_type)
    save_results(sol, vnames, noise_type)
    
    println("\n✅ Successfully completed: $noise_type\n")
end

# Compare all results
compare_noise_types(results, var_names)

println("\n" * "="^80)
println("ALL SIMULATIONS COMPLETE!")
println("="^80)
println("\nResults saved in:")
println("  - imgs/bio/symbolic/")
println("  - results/symbolic/")
toc() 

verify_noise_differences(results)

# ============================================================================
# This will verify that the noise types are actually different
# ============================================================================

println("\n" * "="^80)
println("DETAILED VERIFICATION OF NOISE DIFFERENCES")
println("="^80)

using Statistics

# Function to calculate power spectral density
function simple_psd(data, dt=1.0)
    n = length(data)
    data_centered = data .- mean(data)
    fft_result = abs.(fft(data_centered))
    freqs = fftfreq(n, 1/dt)
    return freqs[1:n÷2], fft_result[1:n÷2]
end

# Compare the three noise types
noise_comparison = Dict()

for nt in noise_list
    if haskey(results, nt)
        sol, vn = results[nt]
        df = DataFrame(sol)
        
        # FIX: Check if we have enough data points
        if nrow(df) < 100
            println("⚠️  Warning: $nt has only $(nrow(df)) points, skipping comparison")
            continue
        end
        
        dsol_window = filter(:timestamp => t -> t ≥ (maximum(df.timestamp) - 1000), df)
        caer_idx = df_find_column(dsol_window, "caer")
        cac_idx = df_find_column(dsol_window, "cac")
        
        if !isnothing(caer_idx) && !isnothing(cac_idx)
            # FIX: Safe indexing - take last available points
            n_points = min(1000, nrow(df))
            start_idx = max(1, nrow(df) - n_points + 1)

            caer_data = df[start_idx:end, caer_idx]
            cac_data = df[start_idx:end, cac_idx]
            
            # Filter out NaN/Inf values
            valid_mask = isfinite.(caer_data) .& isfinite.(cac_data)
            if sum(valid_mask) < 10
                println("⚠️  Warning: $nt has insufficient valid data")
                continue
            end
            
            caer_data = caer_data[valid_mask]
            cac_data = cac_data[valid_mask]
            
            noise_comparison[nt] = Dict(
                :caer_mean => mean(caer_data),
                :caer_std => std(caer_data),
                :cac_mean => mean(cac_data),
                :cac_std => std(cac_data),
                :correlation => cor(caer_data, cac_data),
                :caer_data => caer_data,
                :cac_data => cac_data
            )
            
            println("\n$nt:")
            println("  Valid points: $(length(caer_data))")
            println("  CaER: μ=$(round(mean(caer_data), digits=6)), σ=$(round(std(caer_data), digits=6))")
            println("  CaC:  μ=$(round(mean(cac_data), digits=6)), σ=$(round(std(cac_data), digits=6))")
            println("  CaER-CaC correlation: $(round(cor(caer_data, cac_data), digits=4))")
        end
    end
end

# Check if they're actually different
println("\n" * "-"^80)
println("PAIRWISE COMPARISONS:")
println("-"^80)

all_identical = true

# Use a non-constant local name to avoid clashing with any existing `n`
N = length(noise_list)
corrmat = fill(NaN, N, N)
maemat  = fill(NaN, N, N)

for i in 1:N
    corrmat[i,i] = 1.0
    maemat[i,i]  = 0.0
    for j in (i+1):N
        nt1, nt2 = noise_list[i], noise_list[j]
        if haskey(noise_comparison, nt1) && haskey(noise_comparison, nt2)
            d1 = noise_comparison[nt1][:caer_data]
            d2 = noise_comparison[nt2][:caer_data]

            # align by taking the last min-length segment
            minlen = min(length(d1), length(d2))
            if minlen < 5
                # leave NaN entries for insufficient overlap
                continue
            end
            d1s = d1[end-minlen+1:end]
            d2s = d2[end-minlen+1:end]

            c = cor(d1s, d2s)
            m = mean(abs.(d1s .- d2s))

            corrmat[i,j] = c
            corrmat[j,i] = c
            maemat[i,j]  = m
            maemat[j,i]  = m

            # keep a simple summary flag consistent with previous logic
            if c <= 0.95
                all_identical = false
            end
        end
    end
end

# Nicely print matrices
using Printf

println("\nPairwise Correlation matrix (rows/cols = $(join(string.(noise_list), ", ")))")
# header
print(rpad("", 18))
for t in noise_list
    print(rpad(string(t), 12))
end
println()
# rows
for i in 1:N
    print(rpad(string(noise_list[i]), 18))
    for j in 1:N
        v = corrmat[i,j]
        if isnan(v)
            print(rpad("-", 12))
        else
            print(rpad(@sprintf("%6.4f", v), 12))
        end
    end
    println()
end

println("\nPairwise Mean Absolute Error (MAE) matrix (same ordering):")
print(rpad("", 18))
for t in noise_list
    print(rpad(string(t), 12))
end
println()
for i in 1:N
    print(rpad(string(noise_list[i]), 18))
    for j in 1:N
        v = maemat[i,j]
        if isnan(v)
            print(rpad("-", 12))
        else
            print(rpad(@sprintf("%8.4e", v), 12))
        end
    end
    println()
end

# keep the boolean summary output consistent with the rest of the script
if all_identical
    println("\nStatus: ⚠️  All pairwise correlations are > 0.95 (possible identical/noise issue)")
else
    println("\nStatus: ✅ At least one pair differs (correlation ≤ 0.95)")
end

# Create a visualization comparing all three
println("\nCreating comparison visualization...")

# Build stacked comparison plots per-noise (safe column lookup per DataFrame)
plots = []
for nt in noise_list
    if haskey(results, nt)
        sol, _ = results[nt]
        df = DataFrame(sol)
        window = filter(:timestamp => t -> t ≥ (maximum(df.timestamp) - 200), df)

        # Find caer column inside this window DataFrame (safe for different column orders)
        caer_col = df_find_column(window, "caer")
        if isnothing(caer_col)
            println("Warning: could not find 'caer' column for $nt, skipping")
            continue
        end

        p = plot(window.timestamp, window[!, caer_col],
                 xlabel="Time (s)", ylabel="CaER (μM)",
                 title=string(nt), legend=false, linewidth=1.5)
        push!(plots, p)
    end
end

if !isempty(plots)
    p_compare = plot(plots..., layout=(length(plots), 1), size=(1200, 900))
    mkpath("imgs/bio/symbolic/")
    savefig(p_compare, "imgs/bio/symbolic/noise_type_comparison.png")
    println("Saved: imgs/bio/symbolic/noise_type_comparison.png")
else
    println("No CaER plots created: no valid data for selected noise types.")
end

# Overlay of the same windows (use per-window column lookup)
p_overlay = plot(xlabel="Time (s)", ylabel="CaER (μM)",
                 title="All Noise Types Overlaid", legend=:topright, size=(1200, 600))

ci = 1
for nt in noise_list
    if haskey(results, nt)
        sol, _ = results[nt]
        df = DataFrame(sol)
        window = filter(:timestamp => t -> t ≥ (maximum(df.timestamp) - 200), df)

        caer_col = df_find_column(window, "caer")
        if isnothing(caer_col)
            println("Warning: could not find 'caer' column for $nt, skipping overlay")
            continue
        end

        plot!(p_overlay, window.timestamp, window[!, caer_col],
              label=string(nt), color=colors[ci <= length(colors) ? ci : 1],
              linewidth=2, alpha=0.7)
        ci += 1
    end
end

# Save overlay only if at least one series was plotted
if !isempty(p_overlay.series_list)
    savefig(p_overlay, "imgs/bio/symbolic/noise_overlay_comparison.png")
    println("Saved: imgs/bio/symbolic/noise_overlay_comparison.png")
else
    println("No overlay saved: no valid CaER series found.")
end

println("\n" * "="^80)
println("VERIFICATION COMPLETE")
println("="^80)


if all_identical
    println("\n" * "="^80)
    println("DEBUG: Checking noise matrix structure")
    println("="^80)
    
    for nt in noise_list
        println("\nRebuilding $nt to check noise matrix...")
        all_eqs, noiseeqs, u0 = build_stochastic_model(nt)
        
        println("  Noise matrix size: $(size(noiseeqs))")
        println("  Number of non-zero entries: $(count(!iszero, noiseeqs))")
        
        # Find which equations have noise
        non_zero_eqs = findall(row -> any(!iszero, noiseeqs[row, :]), 1:size(noiseeqs, 1))
        println("  Equations with noise: $(length(non_zero_eqs))")
        
        if length(non_zero_eqs) > 0
            println("  First few non-zero coefficients:")
            for idx in non_zero_eqs[1:min(5, length(non_zero_eqs))]
                println("    Eq $idx: $(noiseeqs[idx, :])")
            end
        end
    end
end

# ============================================================================
# NEW FUNCTION: Scan IP3 and Plot Differences
# ============================================================================

function scan_ip3_experiment()
    println("\n" * "="^80)
    println("RUNNING IP3 SCAN (0.1 - 2.0) FOR ALL NOISE TYPES")
    println("="^80)

    target_noise_types = [:none, :additive, :multiplicative, :state_dependent, :jump]
    ip3_values = 0.1:0.2:2.0
    
    # Store results: Dict[NoiseType -> DataFrame]
    results_all = Dict()
    diagnostics = Dict()

    for noise_type in target_noise_types
        println("\nProcessing $noise_type...")
        
        # Initialize DataFrame with String column names
        df_res = DataFrame()
        df_res[!, "ip3"] = Float64[]
        for v in important_variables
            df_res[!, v] = Float64[]
        end
        
        # Diagnostics storage
        diag = Dict(
            :converged => Int[],
            :final_time => Float64[],
            :n_steps => Int[],
            :solve_time => Float64[]
        )
        
        # Loop IP3
        for val in ip3_values
            print("  IP3 = $val : ")
            
            # Time the solve
            t_start = time()
            
            try
                # Solve with current IP3 value
                sol, _ = simulate_model(noise_type; tspan=(0.0, 4000.0), ip3_val=val)
                
                t_solve = time() - t_start
                
                # Check convergence
                converged = (Symbol(sol.retcode) == :Success || Symbol(sol.retcode) == :Default) ? 1 : 0
                final_t = sol.t[end]
                n_steps = length(sol.t)
                
                push!(diag[:converged], converged)
                push!(diag[:final_time], final_t)
                push!(diag[:n_steps], n_steps)
                push!(diag[:solve_time], t_solve)
                
                println("✓ ($(round(t_solve, digits=2))s, $n_steps steps, status=$(sol.retcode))")
                
                # Convert to DataFrame
                df_sol = DataFrame(sol)
                n_rows = nrow(df_sol)
                
                println("    DEBUG: DataFrame has $n_rows rows, $(ncol(df_sol)) columns")
                println("    DEBUG: Column names: $(names(df_sol)[1:min(5, ncol(df_sol))])...")
                
                # Calculate steady state (last 20% of trajectory)
                if n_rows > 50 && converged == 1  # Increased threshold for safety
                    start_idx = floor(Int, n_rows * 0.8)
                    row_data = Dict{String, Float64}("ip3" => val)
                    
                    n_found = 0
                    n_missing = 0
                    
                    for var in important_variables
                        col_idx = df_find_column(df_sol, var)
                        if !isnothing(col_idx)
                            data_segment = df_sol[start_idx:end, col_idx]
                            valid_data = filter(!isnan, data_segment)
                            if length(valid_data) > 0
                                row_data[var] = mean(valid_data)
                                n_found += 1
                            else
                                row_data[var] = NaN
                                n_missing += 1
                            end
                        else
                            println("    WARNING: Could not find column for '$var'")
                            row_data[var] = NaN
                            n_missing += 1
                        end
                    end
                    
                    println("    ✓ Extracted: $n_found/$n_missing variables found/missing")
                    push!(df_res, row_data)
                    
                else
                    println("    ⚠️  Insufficient data: n_rows=$n_rows, converged=$converged")
                    fallback_row = Dict{String, Float64}("ip3" => val)
                    for v in important_variables
                        fallback_row[v] = NaN
                    end
                    push!(df_res, fallback_row)
                end
                
            catch e
                println("    ❌ ERROR: $e")
                t_solve = time() - t_start
                
                # Record failure
                push!(diag[:converged], 0)
                push!(diag[:final_time], 0.0)
                push!(diag[:n_steps], 0)
                push!(diag[:solve_time], t_solve)
                
                # Add NaN row
                fallback_row = Dict{String, Float64}("ip3" => val)
                for v in important_variables
                    fallback_row[v] = NaN
                end
                push!(df_res, fallback_row)
            end
        end
        
        results_all[noise_type] = df_res
        diagnostics[noise_type] = diag
        
        # Print summary statistics
        println("\n  Summary for $noise_type:")
        println("    Convergence rate: $(sum(diag[:converged]))/$(length(diag[:converged]))")
        if length(diag[:solve_time]) > 0
            println("    Avg solve time: $(round(mean(diag[:solve_time]), digits=2))s")
            println("    Avg steps: $(round(mean(filter(x -> x > 0, diag[:n_steps])), digits=0))")
        end
    end

    # Print detailed diagnostics
    println("\n" * "="^80)
    println("DETAILED DIAGNOSTICS")
    println("="^80)
    
    for (nt, diag) in diagnostics
        println("\n$nt:")
        n_converged = sum(diag[:converged])
        n_total = length(diag[:converged])
        println("  Convergence: $n_converged/$n_total successful")
        
        if n_converged > 0
            successful_times = diag[:solve_time][diag[:converged] .== 1]
            successful_steps = diag[:n_steps][diag[:converged] .== 1]
            
            println("  Solve time: min=$(round(minimum(successful_times), digits=2))s, " *
                    "max=$(round(maximum(successful_times), digits=2))s, " *
                    "mean=$(round(mean(successful_times), digits=2))s")
            println("  Steps: min=$(minimum(successful_steps)), " *
                    "max=$(maximum(successful_steps)), " *
                    "mean=$(round(mean(successful_steps), digits=0))")
            
            # Check stability
            successful_finals = diag[:final_time][diag[:converged] .== 1]
            stable = all(successful_finals .>= 3900.0)
            println("  Stability: $(stable ? "✓ All reached t=4000" : "⚠️  Some terminated early")")
            if !stable
                early_indices = findall((diag[:final_time] .< 3900.0) .& (diag[:converged] .== 1))
                println("    Early terminations at IP3 = $(ip3_values[early_indices])")
            end
        else
            println("  ⚠️  No successful runs to analyze")
        end
    end

    # Plotting Comparison (rest remains the same)
    println("\n" * "="^80)
    println("GENERATING COMPARISON PLOTS")
    println("="^80)
    mkpath("imgs/bio/scan/")
    
    for var in important_variables
        p = plot(xlabel="IP3 (μM)", ylabel="$var (steady state)", 
                 title="Steady State $var vs IP3", 
                 legend=:outertopright, size=(800, 600))
        
        color_map = Dict(
            :none => :black,
            :additive => :blue,
            :multiplicative => :red,
            :state_dependent => :green,
            :jump => :purple
        )
        
        for (nt, df) in results_all
            if nrow(df) > 0 && var in names(df)
                valid_idx = .!isnan.(df[!, var])
                if sum(valid_idx) > 0
                    plot!(p, df[!, "ip3"][valid_idx], df[!, var][valid_idx], 
                          label=string(nt), 
                          marker=:circle, 
                          markersize=5,
                          linewidth=2,
                          color=get(color_map, nt, :auto))
                end
            end
        end
        
        savefig(p, "imgs/bio/scan/compare_$(var)_ip3.png")
        println("Saved: imgs/bio/scan/compare_$(var)_ip3.png")
    end
    
    # Combined plot
    println("\nCreating combined overview plot...")
    plots_combined = []
    for var in important_variables
        p = plot(xlabel="IP3 (μM)", ylabel=var, title=var, legend=false, size=(400, 300))
        for (nt, df) in results_all
            if nrow(df) > 0 && var in names(df)
                valid_idx = .!isnan.(df[!, var])
                if sum(valid_idx) > 0
                    plot!(p, df[!, "ip3"][valid_idx], df[!, var][valid_idx], 
                          linewidth=2, marker=:circle, markersize=3)
                end
            end
        end
        push!(plots_combined, p)
    end
    
    p_all = plot(plots_combined..., layout=(2, 3), size=(1400, 800))
    savefig(p_all, "imgs/bio/scan/all_variables_overview.png")
    println("Saved: imgs/bio/scan/all_variables_overview.png")
    
    println("\n" * "="^80)
    println("IP3 SCAN COMPLETE")
    println("="^80)
    
    return results_all, diagnostics
end

# Run the experiment
scan_results, scan_diagnostics = scan_ip3_experiment()