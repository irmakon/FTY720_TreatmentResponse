libname h "D:\Users\mse0jxh8\Documents\My SAS Files\9.4\N2301";
libname n2301 "E:\SDT-RP-11223\files_NVT_SA_FTY720D2301\Files\Analysis Ready Datasets\SAS_analysis";

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
	          other = "Other"
				  ;
run;

data rel1;
    set n2301.a_rel (keep = stysid1a RLPCFM1C DAY_1N FLAG_REL RELTM1);
	where RLPCFM1C = 1 and FLAG_REL = 1;
	format FLAG_REL;
run;

proc sort data = rel1;
	by stysid1a day_1n;
run;

data rel2;
    set rel1;
	by stysid1a day_1n;
	if first.stysid1a then output;
run;

data demo;
    set n2301.a_dmg (keep = stysid1a TRTREG1C SEX1C RCE5C age_1n_cat);
*	format RCE5C rce_.;
run;

data ae;
    set n2301.a_aeinf (keep = stysid1a TRTSTT1D SOC_ABBR SOC_TXT AEVSER1C AEVSTT1D DYSTT_1N);
	where SOC_TXT ne "" and DYSTT_1N > 0;               *DYSTT_1N missing = General disorder...;
	trt_stdt = input(TRTSTT1D,date10.);
	ae_stdt = input(AEVSTT1D,date10.);
	diff_trt_ae = ae_stdt - trt_stdt;
	format trt_stdt ae_stdt date10.;
run;

data ae_sae1;
    set ae;
	where aevser1c = 1;
run;

proc sort data = ae_sae1;
	by stysid1a dystt_1n;
run;

data ae_sae2;
    set ae_sae1;
	by stysid1a dystt_1n;
    sae_startday = dystt_1n;
	if first.stysid1a then output;
	keep stysid1a aevser1c sae_startday;
run;

data ae_infec1;
    set ae;
	where soc_abbr = "Infec";
run;

proc sort data = ae_infec1;
	by stysid1a dystt_1n;
run;

data ae_infec2;
    set ae_infec1;
	by stysid1a dystt_1n;
    infec_startday = dystt_1n;
	infec = 1;
	if first.stysid1a then output;
	keep stysid1a infec infec_startday;
run;

data ae_neopl1;
    set ae;
	where soc_abbr = "Neopl";
run;

proc sort data = ae_neopl1;
	by stysid1a dystt_1n;
run;

data ae_neopl2;
    set ae_neopl1;
	by stysid1a dystt_1n;
    neopl_startday = dystt_1n;
	neopl = 1;
	if first.stysid1a then output;
	keep stysid1a neopl neopl_startday;
run;

data visit;
    set n2301.a_vis (keep = stysid1a dyvis_1n visnam1a);
run;

proc sort data = visit;
    by stysid1a dyvis_1n;
run;

data last_visit;
    set visit;
    by stysid1a dyvis_1n;
    if last.stysid1a then output;
	rename dyvis_1n = DayLastVisit;
	label dyvis_1n = "Study day of last visit";
run;

data death;
    set n2301.a_cmp (keep = stysid1a dydth_1n dcnrsn1c);
	where dcnrsn1c = 10;
    death = 1;
	format dcnrsn1c dcnr11_.;
run;

proc sort data = death nodupkey;
	by stysid1a;
run;

data discont_ae;
    set n2301.a_cmp (keep = stysid1a dcnrsn1c dylst_1n);
	where dcnrsn1c = 1;
    discont_ae = 1;
	discont_day = dylst_1n;
	format dcnrsn1c dcnr11_.;
	drop dylst_1n;
run;

proc sort data = discont_ae nodupkey;
	by stysid1a;
run;

data edss;
    set n2301.a_edss (keep = stysid1a FCDPG3N1 TCDPG3N1);
	format _all_;
run;

proc sort data = edss nodupkey;
    by _all_;
run;

data mri1;
    set n2301.a_mri (keep = stysid1a DAYNI_1N MRIRSL1A RSLVAL1N TVOLT2CH FRET2ITT);
	where MRIRSL1A = "T2New_R" and RSLVAL1N > 0; 
*          or MRIRSL1A = "T2_Vol" and TVOLT2CH > 0;
	format _all_;
run;

proc sort data = mri1;
    by stysid1a DAYNI_1N;
run;

data mri2;
    set mri1;
    by stysid1a DAYNI_1N;
	t2_startday = DAYNI_1N;
	if first.stysid1a then output;
	keep stysid1a t2_startday;
run;

data edss2;
    set n2301.a_edss (keep = stysid1a FCDPG3N1 TCDPG3N1 DAYNI_1N);
	format _all_;
run;

proc sort data = edss2;
    by stysid1a DAYNI_1N;
run;

data edss3;
    set edss2;
    by stysid1a DAYNI_1N;
	rename DAYNI_1N = DayLastEdss;
	label DAYNI_1N = "DayLastEdss";
	if last.stysid1a then output;
run;

data lvis;
    merge last_visit edss3;
    by stysid1a;
run;

/*
data trt;
    set n2301.a_trt;
	format _all_;
run;

data a_cmp;
    set n2301.a_cmp (keep = stysid1a dcnrsn1c dystt_1n LSTSTR1D);
	where dcnrsn1c = 1;
run;

data a_aeinf;
    set n2301.a_aeinf (keep = stysid1a ACNTAK3N AEVSTT1D AEVEND1D);
    where ACNTAK3N = 2;
	format _all_;
run;

data cmp_ae;
    merge a_cmp a_aeinf;
    by stysid1a;
run;

data ident (keep = stysid1a firsmd1o firsmd1d dth1o instudy1 diff_death_trtstart death dcnrsn1c);
    set n2301.a_ident;
	first_dt = input(firsmd1d,date10.);
	format _all_;
	format firsmd1o dth1o datetime20.;
    if dth1o ne . then do;
        diff_death_trtstart = (dth1o - firsmd1o) / (3600 * 24) + 1;
        death = 1;
    end;
run;

data edss (keep = stysid1a kmcensor);
    set n2301.a_edss;
	format _all_;
	format kmcensor datetime20.;
run;

data ident_edss;
    merge ident edss;
	by stysid1a;
    diff = kmcensor - firsmd1o;
	first_dt = input(firsmd1d,date10.);
    diff2 = kmcensor - first_dt;
	format first_dt date10.;
run;
*/

data outcome;
	retain stysid1a DayLastVisit aevser1c sae_startday outcome_sae censor_sae day_sae
           infec infec_startday outcome_infec censor_infec day_infec
		   neopl neopl_startday outcome_neopl censor_neopl day_neopl
		   RLPCFM1C FLAG_REL day_1n RELTM1 outcome_relaps censor_relaps day_relaps
           dydth_1n death outcome_death censor_death day_death
		   dcnrsn1c discont_ae discont_day dcnrsn1c outcome_discont_ae censor_discont_ae day_discont_ae
           t2_startday outcome_t2new censor_t2new day_t2new
		   FCDPG3N1 TCDPG3N1 DayLastEdss outcome_edss censor_edss day_edss
    ;
    merge demo last_visit rel2 ae_sae2 ae_infec2 ae_neopl2 death discont_ae mri2 edss3;
	by stysid1a;
	*Serious AE;
    if sae_startday < 766 and aevser1c = 1 then do;
        outcome_sae = 1;
		censor_sae = 1;
		day_sae = sae_startday;
	end;
	else do;
        if DayLastVisit > 765 then do;
		    censor_sae = 0;
    		day_sae = 765;
		end;
		else do;
		    censor_sae = 0;
    		day_sae = DayLastVisit;
		end;
        if DayLastVisit > 675 then do;
	        outcome_sae = 0;
		end;
		else do;
	        outcome_sae = .;
		end;
	end;
	*Infection;
    if infec_startday < 766 and infec = 1 then do;
        outcome_infec = 1;
		censor_infec = 1;
		day_infec = infec_startday;
	end;
	else do;
        if DayLastVisit > 765 then do;
		    censor_infec = 0;
    		day_infec = 765;
		end;
		else do;
		    censor_infec = 0;
    		day_infec = DayLastVisit;
		end;
        if DayLastVisit > 675 then do;
	        outcome_infec = 0;
		end;
		else do;
	        outcome_infec = .;
		end;
	end;
	*Neoplasm;
    if neopl_startday < 766 and neopl = 1 then do;
        outcome_neopl = 1;
		censor_neopl = 1;
		day_neopl = neopl_startday;
	end;
	else do;
        if DayLastVisit > 765 then do;
		    censor_neopl = 0;
    		day_neopl = 765;
		end;
		else do;
		    censor_neopl = 0;
    		day_neopl = DayLastVisit;
		end;
        if DayLastVisit > 675 then do;
	        outcome_neopl = 0;
		end;
		else do;
	        outcome_neopl = .;
		end;
	end;
	*Relaps;
    if RELTM1 < 766 and RELTM1 ne . then do;
        outcome_relaps = 1;
		censor_relaps = 1;
		day_relaps = RELTM1;
	end;
	else do;
        if DayLastVisit > 765 then do;
		    censor_relaps = 0;
    		day_relaps = 765;
		end;
		else do;
		    censor_relaps = 0;
    		day_relaps = DayLastVisit;
		end;
        if DayLastVisit > 675 then do;
	        outcome_relaps = 0;
		end;
		else do;
	        outcome_relaps = .;
		end;
	end;
	*Death;
    if death = 1 and dydth_1n < 766 then do;
        outcome_death = 1;
		censor_death = 1;
		day_death = dydth_1n;
	end;
	else do;
        if DayLastVisit > 765 then do;
		    censor_death = 0;
    		day_death = 765;
		end;
		else do;
		    censor_death = 0;
    		day_death = DayLastVisit;
		end;
        if DayLastVisit > 675 then do;
	        outcome_death = 0;
		end;
		else do;
	        outcome_death = .;
		end;
	end;
	*Discontinued due to AE;
    if discont_ae = 1 and discont_day < 766 then do;
        outcome_discont_ae = 1;
		censor_discont_ae = 1;
		day_discont_ae = discont_day;
	end;
	else do;
        if DayLastVisit > 765 then do;
		    censor_discont_ae = 0;
    		day_discont_ae = 765;
		end;
		else do;
		    censor_discont_ae = 0;
    		day_discont_ae = DayLastVisit;
		end;
        if DayLastVisit > 675 then do;
	        outcome_discont_ae = 0;
		end;
		else do;
	        outcome_discont_ae = .;
		end;
	end;
	*T2 new/enlarged;
    if t2_startday > 0 and t2_startday < 766 then do;
        outcome_t2new = 1;
		censor_t2new = 1;
		day_t2new = t2_startday;
	end;
	else do;
        if DayLastVisit > 765 then do;
		    censor_t2new = 0;
    		day_t2new = 765;
		end;
		else do;
		    censor_t2new = 0;
    		day_t2new = DayLastVisit;
		end;
        if DayLastVisit > 675 then do;
	        outcome_t2new = 0;
		end;
		else do;
	        outcome_t2new = .;
		end;
	end;
	*EDSS;
    if FCDPG3N1 = 1 then do;
	    if TCDPG3N1 < 766 then do;
            outcome_edss = 1;
	    	censor_edss = 1;
		    day_edss = TCDPG3N1;
		end;
		else do;
            outcome_edss = 0;
	    	censor_edss = 0;
		    day_edss = 765;
		end;
	end;
    if FCDPG3N1 = 0 then do;
        if DayLastVisit > 765 then do;
		    censor_edss = 0;
    		day_edss = 765;
		end;
		else do;
		    censor_edss = 0;
    		day_edss = DayLastVisit;
		end;
        if DayLastVisit > 675 then do;
	        outcome_edss = 0;
		end;
		else do;
	        outcome_edss = .;
		end;
	end;
	*outcome #7;
	if outcome_sae=1 or outcome_discont_ae=1 or outcome_death=1 then do;
	    outcome_7 = 1;
		censor_7 = 1;
		day_7 = min(day_sae, day_discont_ae, day_death);
	end;
	else do;
	    censor_7 = censor_sae;   *=censor_discont_ae=censor_death;
		day_7 = day_sae;         *=day_discont_ae=day_death;
		if day_7 > 675 then do;
		    outcome_7 = 0;
		end;
		else do;
		    outcome_7 = .;
		end;
	end;
	*outcome #10;
	if outcome_infec=1 or outcome_neopl=1 then do;
	    outcome_10 = 1;
		censor_10 = 1;
		day_10 = min(day_infec, day_neopl);
	end;
	else do;
	    censor_10 = censor_infec;   *=censor_neopl;
		day_10 = day_infec;         *=day_neopl;
		if day_10 > 675 then do;
		    outcome_10 = 0;
		end;
		else do;
		    outcome_10 = .;
		end;
	end;
	*outcome #11;
	if outcome_edss=1 or outcome_relaps=1 or outcome_7=1 then do;
	    outcome_11 = 1;
		censor_11 = 1;
		day_11 = min(day_edss, day_relaps, day_7);
	end;
	else do;
	    censor_11 = censor_edss;   *=censor_relaps=censor_7;
		day_11 = day_edss;         *=day_relaps=day_7;
		if day_11 > 675 then do;
		    outcome_11 = 0;
		end;
		else do;
		    outcome_11 = .;
		end;
	end;
	drop sex1c rce5c age_1n_cat trtreg1c visnam1a;
	drop aevser1c sae_startday 
         infec infec_startday 
		 neopl neopl_startday 
		 RLPCFM1C FLAG_REL day_1n RELTM1 
         dydth_1n death 
		 dcnrsn1c discont_ae discont_day dcnrsn1c 
         t2_startday 
		 FCDPG3N1 TCDPG3N1 DayLastEdss 
    ;
run;



/*
data h.outcome;
    set outcome;
run;

data h.edss;
    set edss2;
run;

data h.cmp_ae;
    set cmp_ae;
run;

proc datasets library=work nolist ;
    delete ae
           ae_infec1
           ae_infec2
           ae_neopl1
           ae_neopl2
           ae_sae1
           ae_sae2
           death
           demo
           discont_ae
           edss
           edss2
           edss3
           last_visit
           lvis
           mri1
           mri2
           outcome
           rel1
           rel2
           visit
           ;
quit;
run;

*/