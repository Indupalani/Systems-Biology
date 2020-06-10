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
    model.lb(isnan(model.lb)) = 0;
    model.ub(isnan(model.ub)) = 1000;
    model_sol = optimizeCbModel(model);
    growth_single_org(i) = model_sol.f;                 % Growth of a organism individually
    i = i+1;
end

% Creating a community of gut commensal microorganisms (only includes the
% beneficial microorganisms) and community containing pathogen 
nameTags_com = {'model1'; 'model2'; 'model3'; 'model4';'model5'};
commensal_community = createMultipleSpeciesModel({models(1);models(2);models(3);models(4);models(5)},nameTags_com);
nameTags_gut = {'model1'; 'model2'; 'model3'; 'model4';'model5';'pathomodel1';'pathomodel2'};
Gut_community = createMultipleSpeciesModel({models(1);models(2);models(3);models(4);models(5);models(6);models(7)},nameTags_gut);
[commensal_community, Solution_commensal, result_c] = communitysimulate(commensal_community,models,nameTags_com,diet_constraint_community);
[Gut_community,Solution_gut, result_g] = communitysimulate(Gut_community,models,nameTags_gut,diet_constraint_community);
gutcommunity = result_g.BM; commensalcommunity =result_c.BM;
commensalcommuntiy(end+1:end+2,1) =0;
X_1 = {'Commensal community with pathogens', 'Commensal community'};
Y_1 = [gutcommunity, commensalcommunity]; 
bar(transpose(Y_1))
set(gca,'XTickLabel',X)
legend(legends)
Max_growth = result_g.GRmax;

% Influence of antibiotic on the growth of community
antibioticinhibit_proteinsynthesis = [10,5,4.5,4,3.5,3,2.5,2,1.5,1.25,1,0.5,0];
[Growth_rates,result_p1,result_p2,result_p12] = antibioticsimulation(Gut_community, antibioticinhibit_proteinsynthesis);

% Plotting the growth of a commmunity at different antibiotic concentration
col = 'b-g-r-b:g:m-c-';
legends = {'B.bifidum';'B.longum';'L.salivarius';'L.reuteri';'F.prasutnzii';'Y.enterolitica';'L.monocytogenes'; 'Community'};
i=0; n =1;
for n = 1:size(Growth_rates,1)
   if mod(n,length(legends)) == 0
        yyaxis right;
        plot(antibioticinhibit_proteinsynthesis,Growth_rates(n,:),'k-');
        legend(legends);
        xlabel('Resistance to antibiotic');
        ylabel('Growth rate(1/h)');
        savefig(strcat('Community model',int2str(n/8)));
        figure((n/8)+1);
        n = n+1;
        i =0;
   else
        if n==1 
            figure(1);
        end
        colval = (2*i)+1;
        plot(antibioticinhibit_proteinsynthesis,Growth_rates(n,:),col((colval:(colval+1))));
        hold on;
        n = n+1;
        i = i+1;
    end
end
  
% Analysis of flux 
[X,Y,Inc_flux, Inc_Exc, Dec_flux,Dec_Exc] = flux_analysis(result_p1,result_p2,result_p12,Gut_community,Max_growth);
Y = (Y./sum(Y,2))*100;
bar(Y);
set(gca,'XTickLabel',X);
title('Flux Variation');
ylabel('% of Reaction')
legend({'Resistance acquired by Y.enterocolitis', 'Resistance acquired by L.monocytogenes','Resistance acquired by both'})
saveFigure('Flux_variation');

% The figures show the flux changing in each community with respect to the
% growth and antibiotic uptake.



%%Functions needed for the above program
%Function to simulate the growth of community models 
function [community,sol,result] = communitysimulate(community,models,nameTags,diet_constraint)
[community .infoCom, community .indCom] = getMultiSpeciesModelId(community , nameTags);
k=1;
 for k = 1:length(nameTags)
     rxnBiomass(k,1) = strcat(nameTags(k),models(k).rxns(find(models(k).c)));
     k =k +1;
 end
rxnBiomassId = findRxnIDs(community , rxnBiomass); 
community.infoCom.spBm = rxnBiomass;
community.indCom.spBm = rxnBiomassId;
community = useDiet(community,diet_constraint);
[sol,result] = SteadyCom(community);
end


% Function to simulate the growth of community at varying concentration of antibiotic uptake 
function [Growth_rates,result_p1,result_p2,result_p12]= antibioticsimulation(Gut_community, proteinsynthesislevel)
j=1;
biomass_p1 =[]; biomass_p2 =[]; biomass_p12 = [];
Community_p1 =[]; Community_p2 =[]; Community_p12 =[];
flux_p1 = []; flux_p2 = []; flux_p12 =[];
for j = 1: length(proteinsynthesislevel)
    gut_community_p1 = changeRxnBounds(Gut_community,'pathomodel1pbiosynthesis',-proteinsynthesislevel(j),'l');
    gut_community_p1 = changeRxnBounds(Gut_community,'pathomodel1pbiosynthesis',proteinsynthesislevel(j),'u');
    gut_community_p2 = changeRxnBounds(Gut_community,'pathomodel2pbiosynthesis',-proteinsynthesislevel(j),'l');
    gut_community_p2 = changeRxnBounds(Gut_community,'pathomodel2pbiosynthesis',proteinsynthesislevel(j),'u');
    gut_community_p12 = changeRxnBounds(gut_community_p1,'pathomodel2pbiosynthesis',-proteinsynthesislevel(j),'l');
    gut_community_p12 = changeRxnBounds(gut_community_p1,'pathomodel2pbiosynthesis',proteinsynthesislevel(j),'u');
    [sol_p1(j),result_p1(j)] = SteadyCom(gut_community_p1);
    [sol_p2(j),result_p2(j)] = SteadyCom(gut_community_p2);
    [sol_p12(j),result_p12(j)] = SteadyCom(gut_community_p12);
    biomass_p1 =[biomass_p1, result_p1(j).vBM];
    biomass_p2 =[biomass_p2, result_p2(j).vBM];
    biomass_p12 =[biomass_p12, result_p12(j).vBM];
    Community_p1 = [Community_p1, result_p1(j).GRmax];
    Community_p2 = [Community_p2, result_p2(j).GRmax];
    Community_p12 = [Community_p12, result_p12(j).GRmax];
    flux_p1 = [flux_p1, result_p1(j).flux];
    flux_p2 = [flux_p2, result_p2(j).flux];
    flux_p12 = [flux_p12, result_p12(j).flux];
    j=j+1;
end
Growth_rates = [biomass_p1; Community_p1; biomass_p2;  Community_p2; biomass_p12; Community_p12];
end


% Function to analyze the variation in flux
function [X,Y,Common_Inc_flux,Common_Inc_Exc,Common_Dec_flux,Common_Dec_Exc] = flux_analysis(community_1,community_2,community_3,gut_community,maxgrowth)
k=1;
for k =1:3
    if k==1 
        community = community_1;
    elseif k ==2
        community = community_2;
    else
        community = community_3;
    end
i =1;j=1; num = []; fluxchangingrxns = []; difference_in_flux=[];
for i = 1:12
    if(abs(community(i+1).GRmax - community(i).GRmax) > 0.01*maxgrowth)
        num(j) = i+1;
        j =j+1;
    end
    i =i+1;
end
fluxdiffrxn_1 = find(abs(community(num(1)).flux - community(1).flux) > 1e-5);
if length(num) > 1
   fluxdiffrxn_2 = find(abs(community(num(2)).flux - community(num(1)).flux) > 1e-5); 
   fluxdiffrxn_3 = find(abs(community(num(3)).flux - community(num(2)).flux) > 1e-5);  
   fluxdiffrxn_4 = find(abs(community(num(4)).flux - community(num(3)).flux) > 1e-5);
   fluxchangingrxns = fluxdiffrxn_1(find(ismember(fluxdiffrxn_1,fluxdiffrxn_2)));
   fluxchangingrxns = fluxchangingrxns(find(ismember(fluxchangingrxns,fluxdiffrxn_3)));
   fluxchangingrxns = fluxchangingrxns(find(ismember(fluxchangingrxns,fluxdiffrxn_4)));
   fluxchangingrxns(:,2) = community(1).flux(fluxchangingrxns(:,1));
   fluxchangingrxns(:,3) = community(num(1)).flux(fluxchangingrxns(:,1));
   fluxchangingrxns(:,4) = community(num(2)).flux(fluxchangingrxns(:,1));
   fluxchangingrxns(:,5) = community(num(3)).flux(fluxchangingrxns(:,1));
   fluxchangingrxns(:,6) = community(num(4)).flux(fluxchangingrxns(:,1));
   difference_in_flux(:,1) = fluxchangingrxns(:,1);
   l =1;
   for l =1 : length(num)
    difference_in_flux(:,l+1) = fluxchangingrxns(:,l+2) - fluxchangingrxns(:,l+1);
    l = l+1;
   end
else
    fluxchangingrxns = fluxdiffrxn_1;
    fluxchangingrxns(:,2) = community(1).flux(fluxchangingrxns(:,1));
    difference_in_flux = fluxchangingrxns;
end
    difference_in_flux(:,end) = all(difference_in_flux > 0, 2);
    increasing_flux_rxn_no = difference_in_flux((find(difference_in_flux(:,end))),1);
    decreasing_flux_rxn_no = difference_in_flux((find(~difference_in_flux(:,end))),1);
    increasing_flux_rxns = gut_community.rxns(increasing_flux_rxn_no);
    decreasing_flux_rxns = gut_community.rxns(decreasing_flux_rxn_no);
    Exc_rxn_increasing = increasing_flux_rxns(strmatch('EX',increasing_flux_rxns));
    Exc_rxn_decreasing = decreasing_flux_rxns(strmatch('EX',decreasing_flux_rxns));
    
    if k == 1
        increasingfluxrxns_1 = increasing_flux_rxns; decreasingfluxrxns_1 = decreasing_flux_rxns;
        increasingexcrxns_1 = Exc_rxn_increasing; decreasingexcrxns_1 = Exc_rxn_decreasing;
    elseif k == 2
        increasingfluxrxns_2 = increasing_flux_rxns; decreasingfluxrxns_2 = decreasing_flux_rxns;
        increasingexcrxns_2 = Exc_rxn_increasing; decreasingexcrxns_2 = Exc_rxn_decreasing;
    else
        increasingfluxrxns_12 = increasing_flux_rxns; decreasingfluxrxns_12 = decreasing_flux_rxns;
        increasingexcrxns_12 = Exc_rxn_increasing; decreasingexcrxns_12 = Exc_rxn_decreasing;
    end
    k=k+1;
end
    Increasing_rxns = [length(increasingfluxrxns_1),length(increasingfluxrxns_2),length(increasingfluxrxns_12)];
    Decreasing_rxns = [length(decreasingfluxrxns_1),length(decreasingfluxrxns_2),length(decreasingfluxrxns_12)];
    Inc_Exc_rxns =  [length(increasingexcrxns_1),length(increasingexcrxns_2),length(increasingexcrxns_12)];
    Dec_Exc_rxns = [length(decreasingexcrxns_1),length(decreasingexcrxns_2),length(decreasingexcrxns_12)];
    X = {'Increasing rxns','Decresing rxns','Inc. Exc. rxns','Dec.Exc. rxn'};
    Y = [Increasing_rxns; Decreasing_rxns; Inc_Exc_rxns; Dec_Exc_rxns];
    Inc_flux = increasingfluxrxns_1(find(ismember(increasingfluxrxns_1,increasingfluxrxns_12)));
    Common_Inc_flux = Inc_flux(find(ismember(Inc_flux,increasingfluxrxns_2)));
    Inc_Exc = increasingexcrxns_1(find(ismember(increasingexcrxns_1,increasingexcrxns_12)));
    Common_Inc_Exc = Inc_Exc(find(ismember(Inc_Exc,increasingexcrxns_2)));
    Dec_flux = decreasingfluxrxns_1(find(ismember(decreasingfluxrxns_1,decreasingfluxrxns_12)));
    Common_Dec_flux = Dec_flux(find(ismember(Dec_flux,decreasingfluxrxns_2)));
    Dec_Exc = decreasingexcrxns_1(find(ismember(decreasingexcrxns_1,decreasingexcrxns_12)));
    Common_Dec_Exc = Dec_Exc(find(ismember(Dec_Exc,decreasingexcrxns_2)));
end    