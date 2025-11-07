## PARAMETERS ##

# Conservation relations        units                                                                   Source/Authors
uMmM = 1000                     # (uM/mM)
hs = 3600                       # (h/s)

# Relations 
amtot = 15                      # = atpm + adpm (mM)                                                    (Voorslujis et al., 2024)
atot = 3                        # = atpc + adpc (mM)                                                    (Voorslujis et al., 2024)
ckint = 1                       # = cit + isoc + akg + scoa + suc + fum + mal + oaa (mM)                (Cortassa et al., 2003)
nadtot = 0.8                    # = nadm + nadhm (mM)                                                   (Voorslujis et al., 2024)

const ip3 = 0.7                 # mM IP3 concentration in cytosol (0.5-1 mM)                            (Voorslujis et al., 2024)
const F26B = 0.05               # mM                                                                    (Mulukutla et al., 2015)
const GLCex= 25.0               # mM Glucose concentration in the extracellular space (5-25 mM) 
const glutm = 0.8               # mM 
const aspm = 5                  # mM
const fadm = 0.015              # mM

const RTF = 26.712338           # mV     = R*T/F, 
#R is the gas cte (8.314 J/(mol*K)), T is the abs temp (310K), F is Faraday’s cte (96485 C/mol)         (Magnus and Keizer 1997)

#ph= -log10[H+]                 Equation for proton concentration from pH
const pHm = 8.0                 #  pH in the mitochondria                                               (Mulukutla et al., 2015)
const pHc = 7.3                 #  pH in the cytosol                                                    (Mulukutla et al., 2015)   
#hc= 6.31e-5                     #mM     Cytosolic proton concentration                                  (Buckler et al., 1990 and Casey et al., 2010)
hc= 1000*10^(-pHc)              #mM     Cytosolic proton concentration          pHc= 7.3                (Buckler et al., 1990 and Casey et al., 2010)
hm = 1000*10^(-pHm)             #mM     Mitochondria proton concentration       pHm= 8.0                (Buckler et al., 1990 and Casey et al., 2010)

## Concentrations:
#const Ccadp = 0.54              #mM     Cytosolic adenine diphosphate nucleotide concentration         (Mulukutla et al., 2015) 
#const Ccatp= 0.31               #mM     Cytosolic adenine triphosphate nucleotide concentration        (Mulukutla et al., 2015) 
#const tnadc= 0.3                #mM     Mitochondrial pyruvate concentration 
const Ccamp= 0.03               #mM     Cytosolic adenine monophosphate nucleotide concentration        (Mulukutla et al., 2015)
const Mg= 0.7                   #mM     Cytosolic magnesium concentration                               (Mulukutla et al., 2015)
const MgADP= 0.46               #mM     Cytosolic MgADP concentration                                   (Mulukutla et al., 2015)
const MgATP= 2.69               #mM     Cytosolic MgATP concentration                                   (Mulukutla et al., 2015)                        
const Ccpi= 2.5                 #mM     Cytosolic phosphate concentration                               (Mulukutla et al., 2015) 
const pim = 20.0                #mM     Concentration of Inorganic phosphate in mitochondrial matrix    (Magnus and Keizer 1997)
const LAC= 2.5                  #mM     Cytosolic lactate concentration                                 (Farina et al., 2024)
const nadc= 0.299               #mM     Cytosolic NAD concentration (NAD:NADH = 300) 700-200 in animal model    (Fridlyand et al., 2010)
const nadhc= 0.001              #mM     Cytosolic NADH concentration                                    (Tseng et al., 2024)
const Cc23p2g= 3.1              #mM     Cytosolic 2,3-bisphosphoglycerate concentration                 (Mulukutla et al., 2015) 
const Ccala= 0.2                #mM     Cytosolic alanine concentration                                 (Mulukutla et al., 2015) 
const Ccgsh= 0.57               #mM     Cytosolic gluthione concentration                               (Mulukutla et al., 2015) 
const Ccg16p = 0.1              #mM     Cytosolic glucose 1,6-bisphosphate concentration                (Mulukutla et al., 2015) 

const coa = 0.02                #mM     Concentration of coenzyme A                                     (Cortassa et al., 2003) 
const rhof1 = 0.23              #mM     Concentration of ATPase pumps                                   (Voorslujis et al., 2024)
const DpH = -0.8                #       Ph difference beteen cytosol and mito matrix= pHc - pHm         (Casey et al., 2010)
const rhores = 1.00             #mM     Respiration-driven H+ pump concentration                        (Magnus and Keizer 1997) 
const fc = 0.01                 #       fraction of free Ca2+ in cytosol                                (MK 1998a, Fall-Keizer 2001)
const co2= 1.2                  #mM     Cytosolic CO2 concentration                                     (Mulukutla et al., 2015)

## GLYCOLISIS ##

#Glucose Transporter - Jglctr   (Mulukutla et al., 2015) 
V_mf_glut = 7.67 /hs            #mMs-1  Limiting rate of glucose transporter - Forward                  (Mulukutla et al., 2015)
V_mr_glut = 0.767 /hs           #mMs-1  Limiting rate of glucose transporter - Reverse                  (Mulukutla et al., 2015)
#const k_glc_glut = 1.5          #mM     Michaelis constant of neuronal glucose transporter              (Mulukutla et al., 2015)
const k_glc_glut = 2.87         #mM     Michaelis constant of neuronal glucose transporter              (Berndt et al., 2015)

#Hexokinase - Jhx               (Mulukutla et al., 2015)
V_mf_hk=6.3e3 /hs               #mMs-1 Limiting rate of Hexokinase                                      (Mulukutla et al., 2014)
V_mr_hk=41.0  /hs               #mMs-1 Limiting rate of Hexokinase                                      (Mulukutla et al., 2014)
const k_i_mgatp_hk= 1.0         #mM                                                                     (Mulukutla et al., 2015)
const k_i_g6p_hk= 2.0e-2        #mM                                                                         =
const k_m_glc_hk= 0.1           #mM                                                                         =
const k_m_mgatp_hk= 1.0         #mM                                                                         =
const k_i_mgadp_hk= 1.0         #mM                                                                         =
const k_m_g6p_hk= 0.47          #mM                                                                         =
const k_i_g16bp_hk= 3.0e-2      #mM                                                                         =
const k_i_23bpg_hk= 4.0         #mM                                                                         not specified in 2015
const k_i_gsh_hk=3.0            #mM                                                                         =

#Hexokinase - Jhexo             (Berndt et al., 2015)
V_max_hex = 9.36                #mMs-1 Limiting rate of Hexokinase                                      (Berndt et al., 2015)
const k_m_gluc_hexo = 0.043     #mM                                                                     (Berndt et al., 2015)
const k_m_atp_hexo= 0.37        #mM                                                                     (Berndt et al., 2015)
const k_i_atp_hexo= 0.074       #mM                                                                     (Berndt et al., 2015)
const k_i_g6p_hexo= 0.1         #mM                                                                     (Berndt et al., 2015)

#Glucose Phosphate Isomerase - Jgpi (Mulukutla et al., 2015)
V_mf_gpi= 2.4e3 /hs             #mMs-1 Limiting rate of GPI - Forward                                   (Mulukutla et al., 2015)
V_mr_gpi= 2.0e3 /hs             #mMs-1 Limiting rate of GPI - Reverse                                   (Mulukutla et al., 2015)
const k_f_gpi= 9.6e-1           #mM                                                                     (Mulukutla et al., 2015)
const k_r_gpi= 1.23e-1          #mM                                                                     (Mulukutla et al., 2015)

# Glucose-6-phosphate isomerase - Jg6pi (Berndt et al., 2015)
V_max_g6pi = 24.4               #mM  Limiting rate of GPI                                               (Berndt et al., 2015)
const k_e_g6pi = 0.5157         #mM                                                                     (Berndt et al., 2015)
const k_m_g6p_g6pi= 0.593       #mM                                                                     (Berndt et al., 2015)
const k_m_f6p_g6pi= 0.095       #mM                                                                     (Berndt et al., 2015)

# 6-Phosphofructo-2-Kinase - Jpfk2 (Mulukutla et al., 2015)
V_f_pfk2= 300.0 /hs             #mMs-1 Limiting rate of Phosphofructo-2-Kinase                          (Mulukutla et al., 2015)
const k_i_atp_pfk2= 0.15        #mM                                                                     (Mulukutla et al., 2015)
const k_m_f6p_pfk2= 0.032       #mM                                                                     (Mulukutla et al., 2015)
const k_m_atp_pfk2= 0.15        #mM                                                                     (Mulukutla et al., 2015)
const k_m_adp_pfk2= 0.062       #mM                                                                         =
const k_eq_pfk2= 16.0           #mM                                                                         =
const k_m_f26p_pfk2 = 0.02      #mM                                                                         =
const k_i_adp_pfk2= 0.23        #mM                                                                         =
const k_i_f26p_pfk2= 0.008      #mM                                                                         =
const k_i_f6p_pfk2= 0.001       #mM                                                                         =
const k_akt_pfk2= 0.5           #mM                                                                         =
const k_i_pep_pfk2= 0.013       #mM                                                                         =

#Fructose-2,6-Bisphosphatase - Jf26b (Mulukutla et al., 2015)
V_m_f26bp= 13.86 /hs            #mMs-1 Limiting rate of Fructose-2,6-Bisphosphatase                     (Mulukutla et al., 2015)
const k_i_f6p_f26b= 25.0e-3     #mM                                                                     (Mulukutla et al., 2015)
const k_m_f26p_f26b= 1.0e-3     #mM                                                                     (Mulukutla et al., 2015)

# Phosphofructokinase I - Jpfki (Berndt et al., 2015)
V_max_pfki = 49.6               #mMs-1 Limiting rate of Phosphofructokinase                             (Berndt et al., 2015)
const k_m_f6p_pfki = 0.111      #mM                                                                     (Berndt et al., 2015)
const K_0 = 0.55                #mM                                                                         =
const k_a_f26p_pfki = 0.0042    #mM                                                                         =
const n_f26p_pfki = 5.5         #mM                                                                         =
const k_m_atp_pfki = 0.04       #mM                                                                         =
const n_pfki=1.8                #mM                                                                         =
const k_i_atp_pfki=1.2          #mM                                                                         =
const k_a_fru26p=0.005          #mM                                                                         =

#Phosphofructokinase - Jpfk     (Mulukutla et al., 2015)
#V_f_pfk= 15.5e2 /hs            #mMs-1 Limiting rate of Phosphofructokinase  - Forward                  (Mulukutla et al., 2014)
#V_r_pfk= 67.8 /hs              #mMs-1 Limiting rate of Phosphofructokinase  - Reverse                  (Mulukutla et al., 2014)
V_f_pfk= 2.63e2 /hs             #mMs-1 Limiting rate of Phosphofructokinase  - Forward                  (Mulukutla et al., 2015)
V_r_pfk= 11.53 /hs              #mMs-1 Limiting rate of Phosphofructokinase  - Reverse                  (Mulukutla et al., 2015)
const Lpfk=2e-3                 #                                                                           =
const k_m_atp_pfk= 0.1          #mM                                                                         =
const k_m_mg_pfk= 0.2           #mM                                                                         =
const k_m_23bpg_pfk= 0.5        #mM                                                                         =
const k_m_f6p_pfk= 6.0e-2       #mM                                                                         =
const k_m_f16p_pfk= 0.65        #mM
const k_m_amp_pfk= 0.3          #mM                                                                         =
const k_m_g16bp_pfk= 0.1        #mM                                                                         =
const k_m_pi_pfk= 30.0          #mM                                                                         =
const k_m_f26p_pfk= 5.5e-3      #mM                                                                         =
const k_m_mgatp_pfk = 6.8e-2    #mM                                                                         =
const k_m_mgadp_pfk= 0.54       #mM                                                                         =

#Aldolase - Jald                (Mulukutla et al., 2015)
#V_mf_ald= 6.75e2  /hs          #mMs-1 Limiting rate of Aldose  - Forward                               (Mulukutla et al., 2014)
#V_mr_ald= 2.32e3  /hs          #mMs-1 Limiting rate of Aldose  - Reverse                               (Mulukutla et al., 2014)
V_mf_ald= 1.33e2  /hs           #mMs-1 Limiting rate of Aldose  - Forward                               (Mulukutla et al., 2015)
V_mr_ald= 4.57e2  /hs           #mMs-1 Limiting rate of Aldose  - Reverse                               (Mulukutla et al., 2015)
const k_m_f16p_ald= 5.0e-2      #mM                                                                         =
const k_m_gap_ald= 0.189        #mM                                                                         =
const k_i_dhap_ald= 1.1e-2      #mM                                                                         =
const k_i_23bpg_ald= 1.5        #mM                                                                         =
const k_m_dhap_ald= 3.5e-2      #mM                                                                         =
const k_i_f16p_ald= 1.98e-2     #mM                                                                         =

#Triose Phosphate Isomerase - Jtpi (Mulukutla et al., 2015)
V_mf_tpi= 5.1e2 /hs             #mMs-1 Limiting rate of Triose Phosphate Isomerase - Forward            (Mulukutla et al., 2014;2015)
V_mr_tpi= 2.76e3 /hs            #mMs-1 Limiting rate of Triose Phosphate Isomerase - Reverse            (Mulukutla et al., 2015)
#V_mr_tpi= 46.1  /hs            #mMs-1 Limiting rate of Triose Phosphate Isomerase - Reverse            (Mulukutla et al., 2014)
const k_f_tpi = 1.62e-1         #mM                                                                     (Mulukutla et al., 2015)
const k_r_tpi= 4.3e-1           #mM                                                                     (Mulukutla et al., 2015)

# Glyceraldehyde 3-phosphate dehydrogenase - Jgapdh (Berndt et al., 2015)
V_max_gapdh = 72000.0           #mMs-1 Limiting rate of gapdh                                           (Berndt et al., 2015)
const k_e_gapdh = 0.0868        #mM                                                                     (Berndt et al., 2015)
const k_m_nad_gapdh = 0.01      #mM                                                                     (Berndt et al., 2015)
const k_m_grap_gapdh = 0.101    #mM                                                                     (Berndt et al., 2015)
const k_m_pic_gapdh = 3.9       #mM                                                                     (Berndt et al., 2015)
const k_m_nadh_gapdh=0.027      #mM                                                                     (Berndt et al., 2015)
const k_m_bpg13_gapdh= 0.0035   #mM

#Glyceraldehyde 3-Phosphate Dehydrogenase - Jgapdh (Mulukutla et al., 2015)
V_mf_gapd= 5.317e3  /hs         #mMs-1 Limiting rate of gapdh - Forward                                 (Mulukutla et al., 2015)  
#V_mf_gapd= 1e14                #This work
V_mr_gapd= 3.919e3  /hs         #mMs-1 Limiting rate of gapdh - Reverse                                 (Mulukutla et al., 2015)  
const k_i_gap_gapd= 1.59e-16    #mM                                                                         =
const k_ip_gap_gapd= 0.031      #mM                                                                         =
const k_i_b13_gapd= 1.52e-18    #mM                                                                         =
const k_m_b13_gapd = 0.00671    #mM                                                                         =
const k_m_nadh_gapd= 0.0033     #mM                                                                         =
const k_m_gap_gapd= 0.095       #mM                                                                         =
const k_m_nad_gapd= 0.045       #mM                                                                         =
const k_i_pi_gapd= 3.16         #mM                                                                         =
const k_i_nad_gapd= 0.045       #mM                                                                         =
const k_i_nadh_gapd= 0.01       #mM                                                                         =
const k_ip_b13_gapd= 0.001      #mM                                                                         =

#Phosphoglycerate Kinase - Jpgk (Mulukutla et al., 2015)
V_mf_pgk= 5.96e4   /hs          #mMs-1 Limiting rate of Phosphoglycerate Kinase - Forward               (Mulukutla et al., 2015)
V_mr_pgk= 2.39e4   /hs          #mMs-1 Limiting rate of Phosphoglycerate Kinase - Reverse               (Mulukutla et al., 2015)
const k_i_mgadp_pgk= 0.08       #mM                                                                         =
const k_m_mgadp_pgk= 0.1        #mM                                                                         =
const k_m_b13_pgk= 0.002        #mM                                                                         =
const k_i_mgatp_pgk= 0.186      #mM                                                                         =  
const k_m_pg3_pgk= 1.1          #mM                                                                         =
const k_i_b13_pgk= 1.6          #mM                                                                         =
const k_i_pg3_pgk= 0.205        #mM                                                                         =

#Phosphoglycerate Mutase - Jpgm (Mulukutla et al., 2015)
V_mf_pgm = 4.894e5  /hs         #mMs-1 Limiting rate of Phosphoglycerate Mutase - Forward               (Mulukutla et al., 2015)
V_mr_pgm = 4.395e5  /hs         #mMs-1 Limiting rate of Phosphoglycerate Mutase - Reverse               (Mulukutla et al., 2015)
const k_m_pg3_pgm= 0.168        #mM                                                                         =
const k_m_pg2_pgm= 0.0256       #mM                                                                         =

#Enolase - Jeno                 (Mulukutla et al., 2015)
V_mf_eno= 2.106e4  /hs          #mMs-1 Limiting rate of Enolase - Forward                               (Mulukutla et al., 2015)
V_mr_eno= 5.542e3  /hs          #mMs-1 Limiting rate of Enolase - Reverse                               (Mulukutla et al., 2015)
const k_i_mg_eno= 0.14          #mM                                                                         =
const k_m_pg2_eno= 0.046        #mM                                                                         =
const k_m_pep_eno= 0.11         #mM                                                                         =   
const k_i_pg2_eno=0.046         #mM                                                                         =
const k_i_pep_eno=0.11          #mM                                                                         =

# Pyruvate kinase - Jpk         (Berndt et al., 2015)   
V_max_pkb = 23.76               #mMs-1 Limiting rate of Pyruvate kinase                                 (Berndt et al., 2015) 
const k_m_pep_pkb = 0.074       #mM                                                                     (Berndt et al., 2015) 
const k_m_adp_pkb= 0.42         #mM                                                                     (Berndt et al., 2015) 
const k_i_atp_pkb = 4.4         #mM                                                                     (Berndt et al., 2015) 

# Lactate Dehydrogenase- Jldh   (Mulukutla et al., 2015)
V_mf_ldh= (8.66e2)/hs           #mMs-1 Limiting rate of Lactate Dehydrogenase - Forward                 (Mulukutla et al., 2015)
V_mr_ldh= (2.17e2)/hs           #mMs-1 Limiting rate of Lactate Dehydrogenase - Reverse                 (Mulukutla et al., 2015)
const k_i_nadh_ldh=5.45e-3      #mM                                                                         =
const k_i_nad_ldh= 0.503        #mM                                                                         =
const k_m_pyr_ldh= 0.137        #mM                                                                         =
const k_m_nad_ldh= 0.107        #mM                                                                         =
const k_m_lac_ldh= 1.07         #mM                                                                         =
const k_m_nadh_ldh= 7.43e-3     #mM                                                                         =
const k_i_lac_ldh= 7.33         #mM                                                                         =
const k_i_pyr_ldh= 0.228        #mM                                                                         =

#Pyruvate –Hydrogen shuttle - Jpyrh
#V_max_pyrh = (6.67e12) /hs     #mMs-1 Limiting rate of Pyruvate shuttle                                (Mulukutla et al., 2014)
V_max_pyrh = (1.0e13) /hs       #mMs-1 Limiting rate of Pyruvate shuttle                                (Mulukutla et al., 2015)

# Pyruvate dehydrogenase complex - Jpdh_fad (Berndt et al., 2015)
V_max_pdhf = 13.1               #mMs-1 Limiting rate of pdh (FAD)                                       (Berndt et al., 2015) 
const A_max_cam_pdhf  = 1.7     #                                                                       (Berndt et al., 2015) 
const k_a_cam_pdhf  = 0.001     #mM                                                                     (Berndt et al., 2015) 
const k_m_pyr_pdhf  = 0.068     #mM                                                                     (Berndt et al., 2015)          
const k_m_FADm_pdhf = 0.00001   #mM                                                                         =
const k_m_NADm_pdhf  = 0.041    #mM                                                                         =
const k_m_coam_pdhf  = 0.0047   #mM                                                                         =
const k_i_accoa_pdhf = 0.0004   #mM                                                                         =
const Em_fadpdh= 297.0          #mV                                                                         =
const Em_nadpdh = 150           #mV Arbitrary                                                               =
#k_e_pdhf= exp(((2*Em_fadpdh)+(2*Em_nadpdh))/(1000.0*RTF))

##TCA##

# Jcs -Citrate synthase (CS)
Vmaxcs = 104.0                  #mMs-1  Limiting rate of CS                                             (Voorslujis et al., 2024)
const kmaccoa = 1.2614e-2       #mM     Michaelis constant of CS for AcCoa                              (Cortassa et al.,2003) Dudycha 2000: 1.26e-2 mM; BRENDA: 5-10 uM, Shepherd-Garland 1969: 16 uM)
const kiaccoa = 3.7068e-2       #mM     Inhibition constant of CS for AcCoA                             (Cortassa et al.,2003 or Dudycha 2000) ?
const ksaccoa = 8.0749e-2       #mM     Binding constant of citrate synthase for AcCoA                  (Dudycha 2000)
const kmoaa = 5.0e-3            #mM     Michaelis constant of CS for oxaloacetate                       (Bernart et al., 2015; Matsuoka and Srere; Kurz et al.)
#Jaco -Aconitase (ACO)
const keaco = 0.067             #       Equilibrium constant of aconitase                               (Berndt et al.,2015; Garret & Grisham) equilibrator pH=8, I=0.12M, pMg=3.4, 
const kfaco = 12.5              #s-1    Forward rate constant of aconitase (1/s)                        (Cortassa et al., 2003) 
#Jidh -Isocitrate dehydrogenase (IDH)
Vmaxidh = 1.7767                #mMs-1      Limiting rate of IDH                                        (Cortassa et al., 2003) push conditions
const kmisoc = 1.52             #mM         Michaelis constant of IDH for isocitrate                    (Cortassa et al., 2003)
const kmnad = 0.923             #mM         Michaelis constant of IDH for NAD                           (Cortassa et al., 2003)
const kinadh = 0.19             #mM         Inhibition constant of IDH for NADHm                        (Cortassa et al., 2003)
const kaadp = 6.2e-2            #mM         Activation constant of IDH for ADPm                         (Dudycha 2000 and Cortassa et al., 2003)
const kaca = 1.41               #uM         Activation constant of IDH for mitochondrial Ca2+           (Cortassa et al., 2003)
const kh1 = 8.1e-5              #mM         First ionization constant of IDH                            (Dudycha 2000 and Cortassa et al., 2003)
const kh2 = 5.98e-5             #mM         Second ionisation constant of IDH                           (Dudycha 2000 and Cortassa et al., 2003)
const ni = 2.0                  #           Hill coefficient of IDH for isocitrate                      (Wei et al., 2011 and Cortassa et al., 2011) Cooperativity of Isoc in IDH (lacking in Cortassa 2003)
#Jkgdh -Alpha-ketoglutarate dehydrogenase (alfa-kg) 
Vmaxkgdh = 2.5                  #mMs-1      Limiting rate of aKG dehydrogenase                          (Cortassa et al., 2003)
const kmakg = 1.94              #mM         Michaelis constant of KGDH for aKG                          (Cortassa et al., 2003)
const kmnadkgdh = 0.0387        #mM         Michaelis constant of KGDH for NAD                          (Voorslujis et al., 2024)
const kdca = 1.27               #uM         Dissociation constant of KGDH for mitochondrial Ca2+        (Cortassa et al., 2003)
const kdmg = 0.0308             #mM         Dissociation constant of KGDH for Mg2+                      (Cortassa et al., 2003)
const nakg = 1.2                #           Hill coefficient of KGDH for aKG                            (Cortassa et al., 2003) 
#Jsl -Succinyl-CoA Lyaseor Succinyl-CoA Synthetase (SL)
const kesl = 0.724              #           Equilibrium constant of succinyl coa lyase                  (Flamholz et al., 2012)  equilibrator pH=8, I=0.12M, pMg=3.4 (reac with ADP/ATP: 0.724, with GDP/GTP: 2.152)
const kfsl = 0.127              #mM-1s-1    Forward rate constant of succinyl coa lyase                 (Cortassa et al., 2003)
#Jsdh -Succinate dehydrogenase (SDH)
Vmaxsdh = 0.5                   #mMs-1      Limiting rate of succinate dehydrogenase                    (Cortassa et al., 2003)
const kmsuc = 3.0e-2            #mM         Michaelis constant of SDH for succinate                     (Cortassa et al., 2003)
const kioaasdh = 0.15           #mM         Inhibition constant of SDH for OAA                          (Cortassa et al., 2003)
const kifum = 1.3               #mM         Inhibition constant of SDH for fumarate                     (Cortassa et al., 2003)
#Jfh -Fumarate hydratase (FH)
const kefh = 3.942              #           Equilibrium constant of FH                                  (Flamholz et al., 2012)  equilibrator pH=8, I=0.12M, pMg=3.4
const kffh = 8.3                #s-1        Forward rate constant of FH                                 (Cortassa et al., 2003)
#Jmdh -Malate dehydrogenase (MDH)
Vmdh = 128.0                    # mMs-1     Limiting rate of MDH                                        (Voorslujis et al., 2024)
const Keqmdh = 2.756e-5         #           Equilibrium constant of MDH                                 (Flamholz et al., 2012) equilibrator pH=8, I=0.12M, pMg=3.4 #other values of Keqmdh: #Berndt 2015 : e4 
const kmmal = 0.145             #mM         Michaelis constant of MDH for malate                        (Berndt et al.,2012) #other values of kmmal: #Berndt 2015 : 0.77 
const kmnadmdh = 0.06           #mM         Michaelis constant of MDH for NAD                           (Berndt et al.,2012) #other values of kmnadmdh: #Berndt 2015 : 0.05
const kmoaamdh = 0.017          #mM         Inhibition constant of MDH for OAA                          (Berndt et al.,2012) #other values of kmoaamdh: #Berndt 2015 : 0.04
const kmnadhmdh = 0.044         #mM         Michaelis constant of MDH for NADH                          (Berndt et al.,2012) #other values of kmnadhmdh: #Berndt 2015 : 0.05
#Jf1 - ATPase 
const psiB = 50.0               #mM         Total phase boundary potentials (mV)                        (Magnus and Keizer 1997) Tables 2 and 3
const p1 = 1.346e-8             #                                                                       (Magnus and Keizer 1997) Table 3
const p2 = 7.739e-7             #                                                                       (Magnus and Keizer 1997) Table 3  
const p3 = 6.65e-15             #                                                                       (Magnus and Keizer 1997) Table 3  
const pa = 1.656e-5             #1/s                                                                    (Magnus and Keizer 1997) Table 3  
const pb = 3.373e-7             #1/s                                                                    (Magnus and Keizer 1997) Table 3  
const pc1 = 9.651e-14           #1/s                                                                    (Magnus and Keizer 1997) Table 3  
const pc2 = 4.845e-19           #1/s                                                                    (Magnus and Keizer 1997) Table 3  
const kf1 = 1.71e6              #mM          Equilibrium constant of ATP hydrolysis                     (Cortassa 2003, 2006, 2011)(Magnus and Keizer 1997) Tables 3
#Jhl - M. Proton leak 
const gH = 1.e-5                # mM/(mV * s)    Ionic conductance of the mitochondrial inner membrane      (Cortassa et al., 2003) # Concentrations needed:  DpH = -0.8
#Jhyd - Hydrolisis 
const khyd = 9.e-2              #mMs-1  Hydrolysis rate of ATPc due to cellular activity                    (Voorslujis et al., 2024) Basal hydrolysis rate of ATPc (1/s)
const katpcH = 1.00             #mM     Michaelis constant for ATPc hydrolysis due to cellular activity     (Wacquier et al., 2016)
#Jo- Oxidation/ Respiration  
const r1 = 2.077e-18            #                                                                       (Magnus and Keizer 1997) Table 2                  
const r2 = 1.728e-9             #                                                                       (Magnus and Keizer 1997) Table 2
const r3 = 1.059e-26            #                                                                       (Magnus and Keizer 1997) Table 2
const ra = 6.394e-10            #1/s                                                                    (Magnus and Keizer 1997) Table 2
const rb = 1.762e-13            #1/s                                                                    (Magnus and Keizer 1997) Table 2
const rc1 = 2.656e-19           #1/s                                                                    (Magnus and Keizer 1997) Table 2
const rc2 = 8.632e-27           #1/s                                                                    Magnus and Keizer 1997) Table 2
const g = 0.85                  #       Fitting factor for voltage                                      (Magnus and Keizer 1997) Table 2
const kres = 1.35e18            #                                                                       (Magnus and Keizer 1997) Table 2
#Jant - Adenine nucleotide transporter (Antiporter)     
const vant = 4.0                #mMs-1           Limiting rate of adenine nucleotide translocator       (Voorslujis et al., 2024) #other values of vant: #Fall-Keizer 2001:900 nmol/(mg * min) => 900 * 1000 / 60 = 15 mM/s
const alphac = 0.11111          #                = 0.05 / 0.45                                          (Cortassa et al., 2003)
const alpham = 0.1388889        #                = 0.05 / (0.45 * 0.8)                                  (Cortassa et al., 2003)
const fant = 0.5                #                fraction effective psi                                 (Magnus and Keizer 1997) Table 5 Cortassa 2003
#Jerout - ER LEAK and IP3Rs 
const vip3 = 30.0               #1/s     15. max release rate of Ca2+ through IP3R                      (Voorslujis et al., 2024)
const vleak = 0.15              #1/s     0.15 rate of Ca2+ release by leakage from ER 
const kcai1 = 1.4               #uM      Inhibition constant 1 of IP3R for Ca2+                         Moeien thesis Table 3.1 : 1.3)
const kcai2 = kcai1             #uM      Inhibition constant 1 of IP3R for Ca2+                         Moeien thesis Table 3.1 : 1.3)
const kcaa = 0.70               #uM      Activation constant of IP3R for Ca2+                           Moien 2017 : 0.9
const kip3 = 1.00               #uM      Activation constant of IP3R for IP3                            (Wacquier et al., 2016)
#Jncx
vncx = 2.0e-3                   #mMs-1     Limiting rate of Na+/Ca2+ exchanger                          (Voorslujis et al., 2024)
const psiStar = 91.0            #mV psi offset for Ca2+ transport                                       (Magnus and Keizer 1997) Table 5
const kna = 9.4                 #mM Km (Na+) for NCX                                                    (Magnus and Keizer 1997) Table 5
const kca = 0.375               #uM Km (Ca2+) for NCX                                                   (Bertram et al., 2006) (Cortassa et al., 2003)
const nac = 10.0                #mM Cytosolic Na+ concentration                                         (Voorslujis et al., 2024)
const nam = 5.0                 #mM Mitochondrial Na+ concentration                                     (Voorslujis et al., 2024)
const n = 3.0                   # nb of Na+ binding to NCX (electroneutral/electrogenic: n=2/3) Limiting rate of Na+/Ca2+ exchanger (Magnus and Keizer 1997) Table 5
const b = 0.5                   # NCX dependence on psi (electroneutral/electrogenic: b=0/0.5)          (Magnus and Keizer 1997) Table 5
#Jserva - SERCA 
vserca = 0.12                   #mMs-1 (0.08 y 0.07)   Limiting rate of SERCA                           Moeien thesis 2017 0.455 but multiplied by Ver/Vc~1/3 0.15
const kserca = 0.35             #uM Km of SERCA pumps for Ca2+                                          (Wacquier et al., 2016) under 0.34 and over 0.72 GLY stops being steady
const katpc = 0.05              #mM Km of Serca for ATPc                                                (Wacquier et al., 2016) : 5.e-5 (mM), Moien 2017: 0.06
#Juni - MCU 
const vuni = 0.300              #mMs-1 Limiting rate uniporter                                          (Cortassa et al., 2003): 0.625 mM/s
const L = 110.                  # Keq for uniporter conformations                                       (Magnus and Keizer 1998a)  Table 1
const ma = 2.8                  # uniporter activation cooperativity                                    (Magnus and Keizer 1997)  Table 5
const ktrans = 19.              #uM Kd for uniporter translocated Ca2+                                  (Magnus and Keizer 1998a) Table 1
const kact = 0.38               #uM Kd for uniporter activating Ca2+                                    (Magnus and Keizer 1997) Table 5
#   ADP/ATP exchanger - Jatpex  (Berndt et al., 2015) 
V_max_atpex= 5.4 * 10.0^(-5)    #                                                                       (Berndt et al., 2015) 
const Spsi= 0.3                 #                                                                       (Berndt et al., 2015) 
# Mitochondrial malate dehydrogenase - Jaatm (Berndt et al., 2015) 
V_max_aatm = 32.0               #                                                                       (Berndt et al., 2015) 
const k_e_aatm = 0.147          #                                                                       (Berndt et al., 2015) 

fe = 0.01                       # fraction of free Ca2+ in ER                                           -- Fall-Keizer 2001
fm = 0.0003                     # fraction of free Ca2+ in mitochondria                                 -- Fall-Keizer 2001
alpha = 0.10                    # Ve / Vc -- volumic ratio between ER and cyt
delta = 0.15                    # Vm / Vc -- volumic ration between mito and cyt
cmito = 1.812e-3                # mitochondrial membrane capacitance (mM/mV)
deltaer=0.1                     # fraction of free ca in er * cytosolic to er volum ratio                Moeien thesis deltaer= 0.09
ctot = 1500.0                   # = cac + alpha*caer/fe + delta*cam/fm (uM)
coef1 = fe * ctot / alpha       # coefficients allowing for calculation of caer
coef2 = fe / (alpha * fc)       # coefficients allowing for calculation of caer
coef3 = delta * fe / (alpha * fm) # coefficients allowing for calculation of caer

VDf1B = exp(3.0 * psiB / RTF)      
Pa1 = (pa*(10^(3.0*DpH)))              
Pac1 =  (Pa1 + (pc1 * VDf1B))            
VDresB = exp(6.0*psiB/RTF)         
Ra1 = ra * (10 ^ (6.0*DpH))             
Rac1 = (Ra1 + rc1 * VDresB)               

#scales
s_glctr = 10.0  
s_hx    = 0.8
s_gpi   = 0.6
s_pfk2  = 1.0
s_f26b  = 1.0
s_pfk   = 0.5 # 1 # 0.5
s_ald   = 1.0
s_tpi   = 1.0
s_gapdh = 1.0
s_pgk   = 2.0 #1.5 #2
s_pgm   = 1.0
s_pk    = 2.5 # 1.5
s_ldh   = 1.0
s_pyrh  = 1
s_eno   = 1.0
s_hyd = 0.6 #0.5
s_ant = 0.55
s_f1  = 1.6
s_pdh = 1
