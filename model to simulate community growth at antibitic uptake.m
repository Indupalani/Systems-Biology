clear;
clc;

% Models
%Input the models in the order from the Models folder
Organisms = {'B.bifidum';'B.longum';'L.salivarius';'L.reuteri';'F.prasutnzii';'Y.enterolitica';'L.monocytogenes'};
for i = 1:length(Organisms)
    disp(Organisms(i))
    model = readCbModel();
    models(i) = model;
    i = i+1;
end
 
%Westerndiet
dietConstraints={'EX_fru(e)','-0.14899','1000';'EX_glc_D(e)','-0.14899','1000';'EX_gal(e)','-0.14899','1000';'EX_man(e)','-0.14899','1000';'EX_mnl(e)','-0.14899','1000';'EX_fuc_L(e)','-0.14899','1000';'EX_glcn(e)','-0.14899','1000';'EX_rmn(e)','-0.14899','1000';'EX_arab_L(e)','-0.17878','1000';'EX_drib(e)','-0.17878','1000';'EX_rib_D(e)','-0.17878','1000';'EX_xyl_D(e)','-0.17878','1000';'EX_oxa(e)','-0.44696','1000';'EX_lcts(e)','-0.074493','1000';'EX_malt(e)','-0.074493','1000';'EX_sucr(e)','-0.074493','1000';'EX_melib(e)','-0.074493','1000';'EX_cellb(e)','-0.074493','1000';'EX_tre(e)','-0.074493','1000';'EX_strch1(e)','-0.25734','1000';'EX_amylopect900(e)','-1.5673e-05','1000';'EX_amylose300(e)','-4.7019e-05','1000';'EX_arabinan101(e)','-0.00016628','1000';'EX_arabinogal(e)','-2.1915e-05','1000';'EX_arabinoxyl(e)','-0.00030665','1000';'EX_bglc(e)','-7.05e-08','1000';'EX_cellul(e)','-2.8212e-05','1000';'EX_dextran40(e)','-0.00017632','1000';'EX_galmannan(e)','-1.4106e-05','1000';'EX_glcmannan(e)','-3.2881e-05','1000';'EX_homogal(e)','-0.00012823','1000';'EX_inulin(e)','-0.00047019','1000';'EX_kestopt(e)','-0.0028212','1000';'EX_levan1000(e)','-1.4106e-05','1000';'EX_lmn30(e)','-0.00047019','1000';'EX_lichn(e)','-8.2976e-05','1000';'EX_pect(e)','-3.3387e-05','1000';'EX_pullulan1200(e)','-1.1755e-05','1000';'EX_raffin(e)','-0.0047019','1000';'EX_rhamnogalurI(e)','-1.4492e-05','1000';'EX_rhamnogalurII(e)','-0.00026699','1000';'EX_starch1200(e)','-1.1755e-05','1000';'EX_xylan(e)','-3.2059e-05','1000';'EX_xyluglc(e)','-1.3146e-05','1000';'EX_arachd(e)','-0.003328','1000';'EX_chsterol(e)','-0.004958','1000';'EX_glyc(e)','-1.7997','1000';'EX_hdca(e)','-0.39637','1000';'EX_hdcea(e)','-0.036517','1000';'EX_lnlc(e)','-0.35911','1000';'EX_lnlnca(e)','-0.017565','1000';'EX_lnlncg(e)','-0.017565','1000';'EX_ocdca(e)','-0.16928','1000';'EX_ocdcea(e)','-0.68144','1000';'EX_octa(e)','-0.012943','1000';'EX_ttdca(e)','-0.068676','1000';'EX_ala_L(e)','-1','1000';'EX_cys_L(e)','-1','1000';'EX_ser_L(e)','-1','1000';'EX_arg_L(e)','-0.15','1000';'EX_his_L(e)','-0.15','1000';'EX_ile_L(e)','-0.15','1000';'EX_leu_L(e)','-0.15','1000';'EX_lys_L(e)','-0.15','1000';'EX_asn_L(e)','-0.225','1000';'EX_asp_L(e)','-0.225','1000';'EX_thr_L(e)','-0.225','1000';'EX_glu_L(e)','-0.18','1000';'EX_met_L(e)','-0.18','1000';'EX_gln_L(e)','-0.18','1000';'EX_pro_L(e)','-0.18','1000';'EX_val_L(e)','-0.18','1000';'EX_phe_L(e)','-1','1000';'EX_tyr_L(e)','-1','1000';'EX_gly(e)','-0.45','1000';'EX_trp_L(e)','-0.08182','1000';'EX_12dgr180(e)','-1','1000';'EX_26dap_M(e)','-1','1000';'EX_2dmmq8(e)','-1','1000';'EX_2obut(e)','-1','1000';'EX_3mop(e)','-1','1000';'EX_4abz(e)','-1','1000';'EX_4hbz(e)','-1','1000';'EX_5aop(e)','-1','1000';'EX_ac(e)','-1','1000';'EX_acald(e)','-1','1000';'EX_acgam(e)','-1','1000';'EX_acmana(e)','-1','1000';'EX_acnam(e)','-1','1000';'EX_ade(e)','-1','1000';'EX_adn(e)','-1','1000';'EX_adocbl(e)','-1','1000';'EX_akg(e)','-1','1000';'EX_ala_D(e)','-1','1000';'EX_amet(e)','-1','1000';'EX_amp(e)','-1','1000';'EX_anth(e)','-1','1000';'EX_arab_D(e)','-1','1000';'EX_avite1(e)','-1','1000';'EX_btn(e)','-1','1000';'EX_ca2(e)','-1','1000';'EX_cbl1(e)','-1','1000';'EX_cgly(e)','-1','1000';'EX_chol(e)','-1','1000';'EX_chor(e)','-1','1000';'EX_cit(e)','-1','1000';'EX_cl(e)','-1','1000';'EX_cobalt2(e)','-1','1000';'EX_csn(e)','-1','1000';'EX_cu2(e)','-1','1000';'EX_cytd(e)','-1','1000';'EX_dad_2(e)','-1','1000';'EX_dcyt(e)','-1','1000';'EX_ddca(e)','-1','1000';'EX_dgsn(e)','-1','1000';'EX_etoh(e)','-1','1000';'EX_fald(e)','-1','1000';'EX_fe2(e)','-1','1000';'EX_fe3(e)','-1','1000';'EX_fe3dcit(e)','-1','1000';'EX_fol(e)','-1','1000';'EX_for(e)','-1','1000';'EX_fum(e)','-1','1000';'EX_gam(e)','-1','1000';'EX_glu_D(e)','-1','1000';'EX_glyc3p(e)','-1','1000';'EX_gsn(e)','-1','1000';'EX_gthox(e)','-1','1000';'EX_gthrd(e)','-1','1000';'EX_gua(e)','-1','1000';'EX_h(e)','-1','1000';'EX_h2(e)','-1','1000';'EX_h2s(e)','-1','1000';'EX_hom_L(e)','-1','1000';'EX_hxan(e)','-1','1000';'EX_indole(e)','-1','1000';'EX_ins(e)','-1','1000';'EX_k(e)','-1','1000';'EX_lac_L(e)','-1','1000';'EX_lanost(e)','-1','1000';'EX_mal_L(e)','-1','1000';'EX_metsox_S_L(e)','-1','1000';'EX_mg2(e)','-1','1000';'EX_mn2(e)','-1','1000';'EX_mobd(e)','-1','1000';'EX_mqn7(e)','-1','1000';'EX_mqn8(e)','-1','1000';'EX_na1(e)','-1','1000';'EX_nac(e)','-1','1000';'EX_ncam(e)','-1','1000';'EX_nmn(e)','-1','1000';'EX_no2(e)','-1','1000';'EX_no2(e)','-1','1000';'EX_no3(e)','-1','1000';'EX_orn(e)','-1','1000';'EX_pheme(e)','-1','1000';'EX_phyQ(e)','-1','1000';'EX_pi(e)','-1','1000';'EX_pime(e)','-1','1000';'EX_pnto_R(e)','-1','1000';'EX_ptrc(e)','-1','1000';'EX_pydam(e)','-1','1000';'EX_pydx(e)','-1','1000';'EX_pydx5p(e)','-1','1000';'EX_pydxn(e)','-1','1000';'EX_q8(e)','-1','1000';'EX_retinol(e)','-1','1000';'EX_ribflv(e)','-1','1000';'EX_sel(e)','-1','1000';'EX_sheme(e)','-1','1000';'EX_so4(e)','-1','1000';'EX_spmd(e)','-1','1000';'EX_succ(e)','-1','1000';'EX_thf(e)','-1','1000';'EX_thm(e)','-1','1000';'EX_thymd(e)','-1','1000';'EX_ura(e)','-1','1000';'EX_uri(e)','-1','1000';'EX_vitd3(e)','-1','1000';'EX_xan(e)','-1','1000';'EX_zn2(e)','-1','1000';'EX_meoh(e)','-10','1000';'EX_o2(e)','-10','1000'};
diet_constraint_community = dietConstraints;                        
diet_constraint_community(:,1) = strrep(diet_constraint_community(:,1),'(e)','[u]'); 

% Growth of single organism
i=1;
for i = 1: length(Organisms)
    model = useDiet(models(i),dietConstraints);         % Applying diet constraints to the model
    model.lb(find(isnan(model.lb))) = 0;
    model.ub(find(isnan(model.ub))) = 1000;
    model_sol = optimizeCbModel(model);
    growth_single_org(i) = model_sol.f;                 % Growth of a organism individually
    i = i+1;
end

% Creating a community of gut commensal microorganisms (only includes the
% beneficial microorganisms)
nameTagsModel_com = {'model1'; 'model2'; 'model3'; 'model4';'model5'};
commensal_community = createMultipleSpeciesModel({models(1);models(2);models(3);models(4);models(5)},nameTagsModel_com);
[commensal_community .infoCom, commensal_community .indCom] = getMultiSpeciesModelId(commensal_community , nameTagsModel_com);
 for k = 1:length(nameTagsModel_com)
     rxnBiomass(k,1) = strcat(nameTagsModel_com(k),models(k).rxns(find(models(k).c)));
     k =k +1;
 end
rxnBiomassId = findRxnIDs(commensal_community , rxnBiomass); 
commensal_community.infoCom.spBm = rxnBiomass;
commensal_community.indCom.spBm = rxnBiomassId;
commensal_community = useDiet(commensal_community,diet_constraint_community);
[sol_c,result_c] = SteadyCom(commensal_community);

% Creating a community of gut commensal microorganisms and pathogens
nameTagsModel = {'model1'; 'model2'; 'model3'; 'model4';'model5';'pathomodel1';'pathomodel2'};
gut_community = createMultipleSpeciesModel({models(1);models(2);models(3);models(4);models(5);models(6);models(7)},nameTagsModel);
[gut_community.infoCom, gut_community.indCom] = getMultiSpeciesModelId(gut_community, nameTagsModel);
 for kc = 1:length(Organisms)
     rxnBiomass(kc,1) = strcat(nameTagsModel(kc),models(kc).rxns(find(models(kc).c)));
     kc = kc +1;
 end
rxnBiomassId = findRxnIDs(gut_community, rxnBiomass); 
gut_community.infoCom.spBm = rxnBiomass;
gut_community.indCom.spBm = rxnBiomassId;
gut_community = useDiet(gut_community,diet_constraint_community);
[sol_g,result_g] = SteadyCom(gut_community);

% Influence of antibiotic on the growth of community
antibiotic_uptake = [10,5,4.5,4,3.5,3,2.5,2,1.5,1.25,1,0.5,0];
biomass_p1 =[]; biomass_p2 =[]; biomass_p12 = [];
Community_p1 =[]; Community_p2 =[]; Community_p12 =[];
for j = 1: length(antibiotic_uptake)
    gut_community_p1 = changeRxnBounds(gut_community,'pathomodel1pbiosynthesis',-antibiotic_uptake(j),'l');
    gut_community_p1 = changeRxnBounds(gut_community,'pathomodel1pbiosynthesis',antibiotic_uptake(j),'u');
    gut_community_p2 = changeRxnBounds(gut_community,'pathomodel2pbiosynthesis',-antibiotic_uptake(j),'l');
    gut_community_p2 = changeRxnBounds(gut_community,'pathomodel2pbiosynthesis',antibiotic_uptake(j),'u');
    gut_community_p12 = changeRxnBounds(gut_community_p1,'pathomodel2pbiosynthesis',-antibiotic_uptake(j),'l');
    gut_community_p12 = changeRxnBounds(gut_community_p1,'pathomodel2pbiosynthesis',antibiotic_uptake(j),'u');
    [sol_p1(j),result_p1(j)] = SteadyCom(gut_community_p1);
    [sol_p2(j),result_p2(j)] = SteadyCom(gut_community_p2);
    [sol_p12(j),result_p12(j)] = SteadyCom(gut_community_p12);
    biomass_p1 =[biomass_p1, result_p1(j).vBM];
    biomass_p2 =[biomass_p2, result_p2(j).vBM];
    biomass_p12 =[biomass_p12, result_p12(j).vBM];
    Community_p1 = [Community_p1, result_p1(j).GRmax];
    Community_p2 = [Community_p2, result_p2(j).GRmax];
    Community_p12 = [Community_p12, result_p12(j).GRmax];
    j=j+1;
end   

% Plotting the growth of a commmunity at different antibiotic concentration
col = 'b-g-r-b:g:m-c-';
biomass = [biomass_p1; Community_p1;biomass_p2;  Community_p2;biomass_p12;Community_p12];
legends = {'B.bifidum';'B.longum';'L.salivarius';'L.reuteri';'F.prasutnzii';'Y.enterolitica';'L.monocytogenes'; 'Community'};
figure(1);
i=0; n =1;
for n = 1:size(biomass,1)
   if mod(n, length(legends)) == 0
        yyaxis right;
        plot(antibiotic_uptake,biomass(n,:),'k-');
        legend(legends);
        xlabel('1/antibiotic uptake');
        ylabel('Growth rate(1/h)');
        savefig(strcat('Community model',int2str(n/8)));
        figure((n/8)+1);
        n = n+1;
        i =0;
    else
        colval = (2*i)+1;
        plot(antibiotic_uptake,biomass(n,:),col((colval:(colval+1))));
        hold on;
        n = n+1;
        i = i+1;
    end
end
  

