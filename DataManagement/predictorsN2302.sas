libname h "D:\Users\mse0jxh8\Documents\My SAS Files\9.4\N2302";
libname n2302 "E:\SDT-RP-11223\files_NVT_SA_FTY720D2302\Files\Analysis Ready Datasets\SAS_analysis";

proc format;
    value yno11_ 0 = "no"
	             1 = "yes"
				 ;
    value crn11_ 0 = "no"
	             1 = "yes" 
				 ; *mark if serious;
    value sex11_ 1 = "male"
	             2 = "female"
				 ;
	value dcnr11_ 1 = "Adverse event(s)"
                  2 = "Abnormal laboratory value(s)"
                  3 = "Abnormal test procedure result(s)"
                  4 = "Unsatisfactory therapeutic effect"
                  5 = "Subject’s condition no longer requires study drug"
                  6 = "Protocol violation"
                  7 = "Subject withdrew consent"
                  8 = "Lost to follow-up"
                  9 = "Administrative problems"
                 10 = "Death"
                  ; 
	value rce51_  1 = "Caucasian"
	              2 = "Black"
	              3 = "Asian"
	              7 = "Native American"
	              8 = "Pacific Islander"
	             88 = "Other"
				  ;
	value rce_    1 = "Caucasian"
	          other = "nonCaucasian"
				  ;
    value trtn    1 = "FTY720 1.25mg"
	              2 = "FTY720 0.5mg"
				  3 = "Placebo"
				  ;
run;

data demo;
    set n2302.a_dmg (keep = stysid1a TRTREG1C SEX1C BMI_1N RCE5C age_1n_cat);
	rename sex1c = Sex bmi_1n = BMI rce5c = Race age_1n_cat = Age;
	format RCE5C rce_.;
	label age_1n_cat = "Age category";
run;

data ident;
    set n2302.a_ident (keep = stysid1a trtn priorms TGPDSC1A);
	rename priorms = PriorMSDr trtn = Drug;
run;

data edss;
    set n2302.a_edss (keep = stysid1a vis1n EXDSC1NB FUNBOW1N FUNBRA1N FUNCER1N FUNMEN1N FUNPYR1N
                             FUNSEN1N FUNVSU1N);
	where vis1n = 2; *baseline;
    rename EXDSC1NB = EDSS FUNBOW1N = EDSS_Bowel FUNBRA1N = EDSS_Brain FUNCER1N = EDSS_Cerebellar 
           FUNMEN1N = EDSS_Cerebral FUNPYR1N = EDSS_Pyramidal 
           FUNSEN1N = EDSS_Sensory FUNVSU1N = EDSS_Visual; 
run;

data mri;
    set n2302.a_mri (keep = vis1n stysid1a NUMGD1BA TVOLGDBA TVOLT2BA TVT1HIBA);
	where vis1n = 1; *screening;
	rename NUMGD1BA = MRI_NumGdT1 TVOLGDBA = MRI_VolGdT1 TVOLT2BA = MRI_VolT2 TVT1HIBA = MRI_VolHypT1;
	format _all_;
run;

proc sort data = mri NODUPKEY;
    by _all_;
run;

data mshis;
    set n2302.a_mshis (keep = stysid1a vis1n DURMS1 NUMRLP2N RELTIME);
	where vis1n = 1 and NUMRLP2N ne .; *screening;
	rename DURMS1 = DisDur NUMRLP2N = RlpsNo_2y RELTIME = RlpsDist;
	format _all_;
run;

proc sort data = mshis NODUPKEY;
    by _all_;
run;

data mshis2;
    set n2302.a_mshis (keep = STYSID1A vis1n prvtmt8c prvtmt3c);
	where prvtmt8c ne .;
	format _all_;
run;

proc sort data = mshis2 out = mshis3 NODUPKEY;
    by STYSID1A prvtmt8c;
run;

proc transpose data = mshis3 out = mshis_t;
    by stysid1a;
	id prvtmt8c;
	var prvtmt3c;
run;

data mshis2_t;
    set mshis_t;
	if _1 = 1 or _2 = 1 or _3 = 1 then do;
	    PriorMSDr_IB = 1;
	end;
	else do;
	    PriorMSDr_IB = 0;
	end;
	PriorMSDr_GA = _4;
	PriorMSDr_Nat = _5;
	PriorMSDr_Other = _88;
    PriorMSDr_Count = _1 + _2 + _3 + _4 + _5;
	label PriorMSDr_IB = "Prior Interferon beta-1 treatment"
	      PriorMSDr_GA = "Prior Glatiramer acetate treatment"
	      PriorMSDr_Nat = "Prior Natalizumab treatment"
	      PriorMSDr_Other = "Prior Other MS treatment"
		  PriorMSDr_Count = "Number of prior MS treatments"
		  ;
	drop _NAME_ _LABEL_ _1 _2 _3 _4 _5 _88;
run;
/*
data lrs;
    set n2302.a_lrs (keep = stysid1a vis1n blflg_1c CV_BL_1N PARNAM1C CVUNT_1C LABCAT1C);
	where vis1n = 2 and blflg_1c ="B" and
	      PARNAM1C in ("ABSBAS","ABSEOS","ABSLYM","ABSMON","ABSNEU","ALB","ALKPHS","CREA","DBIL",
                       "GGT","MCH","MCV","PT","SGOT","SGPT","TBIL","UPROTST","WBC");
	format _all_;
run;

data lrs99;
    set n2302.a_lrs (keep = stysid1a CV_BL_1N PARNAM1C);
	where CV_BL_1N ne . and 
	      PARNAM1C in ("ABSBAS","ABSEOS","ABSLYM","ABSMON","ABSNEU","ALB","ALKPHS","CREA","DBIL",
                       "GGT","MCH","MCV","PT","SGOT","SGPT","TBIL","UPROTST","WBC");
	format _all_;
run;

proc sort data = lrs99 nodupkey;
    by stysid1a PARNAM1C CV_BL_1N;
run;

proc means data = lrs99 noprint;
    by stysid1a PARNAM1C;
	var CV_BL_1N;
	output out = out99 n(CV_BL_1N) = nl;
run;                        *-> Baseline-Laborwert gleich pro Patient und Param fuer alle Visiten;

data unique;
    set n2302.a_lrs (keep = PARNAM1C CVUNT_1C LABCAT1C);
	format _all_;
run;

proc sort data = unique nodupkey;
    by _all_;                                                  *-> unique;
run;
*/
data lrs;
    set n2302.a_lrs (keep = stysid1a CV_BL_1N PARNAM1C CVUNT_1C LABCAT1C);
	where CV_BL_1N ne . and 
	      PARNAM1C in ("ABSBAS","ABSEOS","ABSLYM","ABSMON","ABSNEU","ALB","ALKPHS","CREA","DBIL",
                       "GGT","MCH","MCV","SGOT","SGPT","TBIL","UPROTST","WBC");
	format _all_;
run;

proc sort data = lrs nodupkey;
    by _all_;
run;

proc transpose data = lrs out = lrs_t; *Converted Baseline value;
    by stysid1a;
	id PARNAM1C;
	var CV_BL_1N;
run;

%macro renamelab(ds);
    data &ds.2;
    set &ds.;
    %let lab_list = ABSBAS ABSEOS ABSLYM ABSMON ABSNEU ALB ALKPHS CREA DBIL GGT MCH MCV SGOT SGPT TBIL UPROTST WBC
                    ;
	%let b=1;
	%do %while(%scan(&lab_list.,&b.) ne );
	    %let y = %scan(&lab_list.,&b.);
        rename &y. = Lab_&y.;
	    %let b = %eval(&b. + 1);
	%end;
	run;
%mend renamelab;

%renamelab(lrs_t);

/*
data lrs_t2;
    set lrs_t;
	rename ABSBAS = Lab_ABSBAS 
           ABSEOS = Lab_ABSEOS
           ABSLYM = Lab_ABSLYM
           ABSMON = Lab_ABSMON
           ABSNEU = Lab_ABSNEU
           ALB = Lab_ALB
           ALKPHS = Lab_ALKPHS
           CREA = Lab_CREA
           DBIL= Lab_DBIL
           GGT = Lab_GGT
           MCH = Lab_MCH
           MCV = Lab_MCV
           PT = Lab_PT
           SGOT = Lab_SGOT
           SGPT = Lab_SGPT
           TBIL = Lab_TBIL
           UPROTST= Lab_UPROTST
           WBC = Lab_WBC
		   ;
run;
*/

proc transpose data = lrs out = lrs_ct;  *Laboratory test category;
    by stysid1a;
	id PARNAM1C;
	var LABCAT1C;
run;

data lrs_ct2;
    set lrs_ct;
	pos = 1;
	drop stysid1a _name_ _label_;
run;

proc sort data = lrs_ct2;
    by _all_;
run;

data lrs_ct3;
	length ABSBAS ABSEOS ABSLYM ABSMON ABSNEU ALB ALKPHS CREA DBIL GGT MCH MCV SGOT SGPT TBIL UPROTST WBC $10.;
    set lrs_ct2 end = last;
	if last then output;
run;

proc transpose data = lrs out = lrs_ut;  *Standard unit;
    by stysid1a;
	id PARNAM1C;
	var CVUNT_1C;
run;

data lrs_ut2;
    set lrs_ut;
	pos = 2;
	drop stysid1a _name_ _label_;
run;

proc sort data = lrs_ut2;
    by _all_;
run;

data lrs_ut3;
    set lrs_ut2 end = last;
	if last then output;
run;

data lrs_l;  *category (pos=1) + unit (pos=2);
    set lrs_ct3 lrs_ut3;
run;

%macro setlabel(ds,var);  *Label aus category + unit setzen;
    data _null_;
        length cu $20.;
        retain cu;
        set lrs_l;
	    if pos = 1 then do;
	        cu = &var.;
	    end;
	    if pos = 2 then do;
	        cu = catx(' ',cu,&var.);
	    end;
        call symput ('label',cu);
    run;
    data &ds.;
	    set &ds.;
		label lab_&var. = "&label.";
    run;
%mend setlabel;

%macro setlabels(ds);
    %let lab_list = ABSBAS ABSEOS ABSLYM ABSMON ABSNEU ALB ALKPHS CREA DBIL GGT MCH MCV SGOT SGPT TBIL UPROTST WBC
                    ;
	%let b=1;
	%do %while(%scan(&lab_list.,&b.) ne );
	    %let y = %scan(&lab_list.,&b.);
        %setlabel(&ds.,&y.);
	    %let b = %eval(&b. + 1);
	%end;
%mend setlabels;

%setlabels(lrs_t2);

data labor;
    set lrs_t2;
	drop _name_ _label_;
run;

data qol;
    set n2302.a_qol (keep = stysid1a vis1n pag1a ACT4C ANX1C MOB1C PAI4C SFC1C VASBL);
	where vis1n = 2 and pag1a ="EQD_2"; *baseline;
	rename ACT4C = QoL_UsuAct ANX1C = QoL_AnxDep MOB1C = QoL_Mob PAI4C = QoL_Pain
           SFC1C = QoL_Selfcare VASBL = QoL_Vas;
	format _all_;
run;

data msfc;
    set n2302.a_msfc (keep = stysid1a VIS1N HPTBL PASATBL T25FWBL);
	where vis1n = 2; *baseline;
	rename HPTBL = MSFC_9HPT PASATBL = MSFC_PASAT T25FWBL = MSFC_25FW;
	format _all_;
run;

proc sort data = msfc NODUPKEY;
    by _all_;
run;
/*
data cnd99;
    set n2302.a_cnd (keep = stysid1a SOC_ABBR ACTPRB1C SOC_TXT); * SOC_CODE = "" fuer alle;
	where (ACTPRB1C = 0 or ACTPRB1C = .) and SOC_ABBR ne "";
	format _all_;
run;

proc sort data = cnd99 out = cnd99 NODUPKEY;
    by stysid1a SOC_ABBR ACTPRB1C;
run;

proc means data = cnd99 noprint;
    by stysid1a SOC_ABBR;
	var ACTPRB1C;
	output out = out99 n(ACTPRB1C) = nl;
run;            *-> es gibt SOC_ABBR mit ACTPRB1C = 0 und missing pro Patient (und 1);
*/
data cnd;
    set n2302.a_cnd (keep = stysid1a SOC_ABBR ACTPRB1C SOC_TXT); * SOC_CODE = "" fuer alle;
	where SOC_ABBR ne "";
	format _all_;
run;

proc sort data = cnd out = cnd2 NODUPKEY;
    by stysid1a SOC_ABBR ACTPRB1C;
run;

data cnd3;
    set cnd2;
    by stysid1a SOC_ABBR;
    if last.SOC_ABBR then output;
run;

data demo1;
    set demo (keep = stysid1a);
run;

data cnds;
    set cnd (keep = SOC_ABBR);
    stysid1a = "";
run;

proc sort data = cnds nodupkey;
    by SOC_ABBR;
run;

proc transpose data = cnds out = cnds_t;
    id SOC_ABBR;
	var stysid1a;
run;

data cnds_t;
    set cnds_t;
	drop _NAME_;
run;

%macro varlist(ds);
    %local dsid rc nvars i varlist;
	%let dsid = %sysfunc(open(&ds.,is));
	%let nvars = %sysfunc(attrn(&dsid.,nvars));
	%let varlist = ;
	%do i = 1 %to &nvars.;
	    %if %length(&varlist.) eq 0 %then %let varlist = %sysfunc(varname(&dsid.,&i.));
		%else %let varlist = &varlist. %sysfunc(varname(&dsid.,&i.));
	%end;
	%let rc = %sysfunc(close(&dsid.));
	&varlist.
%mend varlist;

%macro demo_cnd;
    %let cnd_list = %varlist(cnds_t);
    data demo2;
        set demo1;
	    %let b=1;
	    %do %while(%scan(&cnd_list.,&b.) ne );
	        %let y = %scan(&cnd_list.,&b.);
            SOC_ABBR = tranwrd("&y.",'_','&');
		    output;
	        %let b = %eval(&b. + 1);
	    %end;
    run;
%mend demo_cnd;

%demo_cnd;

data cnd4;
    merge demo2 (in=in_demo) cnd3 (in=in_cnd);
	by stysid1a SOC_ABBR;
	if in_demo = 1 then demo=1;
	if in_cnd = 1 then cnd=1;
	if in_cnd = 0 and ACTPRB1C = . then ACTPRB1C = 0;
run;

proc transpose data = cnd4 out = cnd_t;
    by stysid1a;
	id SOC_ABBR;
	var ACTPRB1C;
run;

data cnd_t;
    set cnd_t;
	drop _NAME_ _LABEL_;
run;

%macro renamecnd(ds);
    data &ds.2;
    set &ds.;
    %let cnd_list = %varlist(&ds.);
  
	%let b=2;
	%do %while(%scan(&cnd_list.,&b.) ne );
	    %let y = %scan(&cnd_list.,&b.);
        rename &y. = CoDis_&y.;
	    %let b = %eval(&b. + 1);
	%end;
	run;
%mend renamecnd;

%renamecnd(cnd_t);

data cnd_l; *label;
    set cnd;
	keep SOC_ABBR SOC_TXT;
run;

proc sort data = cnd_l nodupkey;
    by _all_;
run;

proc transpose data = cnd_l out = cnd_lt;  *;
	id SOC_ABBR;
	var SOC_TXT;
run;

data cnd_lt;
    set cnd_lt;
	drop _name_ _label_;
run;

%macro setlabel2(ds,var);  *Label fuer CND setzen;
    data _null_;
        set cnd_lt;
        call symput ('label',&var.);
    run;
    data &ds.;
	    set &ds.;
		label CoDis_&var. = "&label.";
    run;
%mend setlabel2;

%macro setlabels2(ds);
    %let cnd_list = %varlist(cnd_lt);
    
	%let b=1;
	%do %while(%scan(&cnd_list.,&b.) ne );
	    %let y = %scan(&cnd_list.,&b.);
        %setlabel2(&ds.,&y.);
	    %let b = %eval(&b. + 1);
	%end;
%mend setlabels2;

%setlabels2(cnd_t2);

data CoDis;
    set cnd_t2;
run;

data cm;
    set n2302.a_cmdatc (keep = stysid1a CONCTIM ATCCODE);
	where ATCCODE ne "" and CONCTIM = 1;
	atc = substr(ATCCODE,1,1);
	format _all_;
run;

proc sort data = cm out = cm2 NODUPKEY;
    by stysid1a atc CONCTIM;
run;

data atccode;
    A = "Alimentary tract and metabolism";
	B = "Blood and blood forming organs";
    C = "Cardiovascular system";
	D = "Dermatologicals";
	G = "Genito urinary system and sex hormones";
	H = "Systemic hormonal preparations, excluding sex hormones and insulins";
	J = "Antiinfective for systemic use";
	L = "Antineoplastic and immunomodulating agents";
	M = "Musculo-skeletal system";
	N = "Nervous system";
	P = "Antiparasitic products, insecticides and repellents";
	R = "Respiratory system";
	S = "Sensory organs";
	V = "Various";
run;

data cm3;
    set cm2;
    by stysid1a atc CONCTIM;
	if last.atc then output;       *ueberfluessig wenn nur CONCTIM = 1;
	keep stysid1a atc CONCTIM;
run;

%macro demo_cm;
    %let cm_list = %varlist(atccode);
    data demo3;
        set demo1;
		length atc $7.;
	    %let b=1;
	    %do %while(%scan(&cm_list.,&b.) ne );
	        %let y = %scan(&cm_list.,&b.);
            atc = "&y.";
		    output;
	        %let b = %eval(&b. + 1);
	    %end;
    run;
%mend demo_cm;

%demo_cm;

data cm4;
    merge demo3 (in=in_demo) cm3 (in=in_cm);
	by stysid1a atc;
	if in_demo = 1 then demo = 1;
	if in_cm = 1 then cm = 1;
	if in_cm = 0 and CONCTIM = . then CONCTIM = 0;
run;

proc transpose data = cm4 out = cm_t;
    by stysid1a;
	id atc;
	var CONCTIM;
run;

%macro renamecm(ds);
    data &ds.2;
    set &ds.;
    %let cm_list = %varlist(atccode);
  
	%let b=1;
	%do %while(%scan(&cm_list.,&b.) ne );
	    %let y = %scan(&cm_list.,&b.);
        rename &y. = CoMed_&y.;
	    %let b = %eval(&b. + 1);
	%end;
	run;
%mend renamecm;

%renamecm(cm_t);

%macro setlabel3(ds,var);  *Label fuer CM setzen;
    data _null_;
        set atccode;
        call symput ('label',&var.);
    run;
    data &ds.;
	    set &ds.;
		label CoMed_&var. = "&label.";
    run;
%mend setlabel3;

%macro setlabels3(ds);
    %let cm_list = %varlist(atccode);
    
	%let b=1;
	%do %while(%scan(&cm_list.,&b.) ne );
	    %let y = %scan(&cm_list.,&b.);
        %setlabel3(&ds.,&y.);
	    %let b = %eval(&b. + 1);
	%end;
%mend setlabels3;

%setlabels3(cm_t2);

data CoMed;
    set cm_t2;
	drop _NAME_ _LABEL_;
run;

data predic;
    merge demo ident edss mri mshis mshis2_t qol msfc;
	by stysid1a;
	drop pag1a vis1n;
run;

/*
data h.labor;
    set labor;
run;

data h.CoDis;
    set CoDis;
run;

data h.predic;
    set predic;
run;

data h.CoMed;
    set CoMed;
run;

data h.atccode;
    set atccode;
run;


proc datasets library=work nolist ;
    delete atccode
           cm
           cm2
           cm3
           cm4
           cm_t
           cm_t2
           comed
           cnd
           cnd2
           cnd3
           cnd4
           cnds
           cnds_t
           cnd_l
           cnd_lt
           cnd_t
           cnd_t2
           codis
           demo
           demo1
           demo2
           demo3
           edss
           ident
           labor
           lrs
           lrs_ct
           lrs_ct2
           lrs_ct3
           lrs_l
           lrs_t
           lrs_t2
           lrs_ut
           lrs_ut2
           lrs_ut3
           mri
           msfc
           mshis
           mshis2
           mshis3
           mshis_t
           mshis2_t
           predic
           qol
           ;
quit;
run;

*/