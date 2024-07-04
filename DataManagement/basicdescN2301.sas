/* Macro zur Erstellung der Basic Description
nummerische Variable: PROC MEANS
Kategorie-/Text-Variable: PROC FREQ
*/
libname n1201 "E:\SDT-RP-11223\files_NVT_SA_FTY720D1201\Files\Analysis Ready Datasets\SAS_analysis";
libname n1201e1 "E:\SDT-RP-11223\files_NVT_SA_FTY720D1201E1\Files\Analysis Ready Datasets\SAS_analysis";
libname n2301 "E:\SDT-RP-11223\files_NVT_SA_FTY720D2301\Files\Analysis Ready Datasets\SAS_analysis";
libname n2302 "E:\SDT-RP-11223\files_NVT_SA_FTY720D2302\Files\Analysis Ready Datasets\SAS_analysis";
libname n2309 "E:\SDT-RP-11223\files_NVT_SA_FTY720D2309\Files\Analysis Ready Datasets\SAS_analysis";
libname r92p1 "E:\SDT-RP-11223\files_RCH-WA21092-PRIM1\Files\Analysis Ready Datasets\SAS_analysis";
libname r92p2 "E:\SDT-RP-11223\files_RCH-WA21092-PRIM2\Files\Analysis Ready Datasets\SAS_analysis";
libname r93p1 "E:\SDT-RP-11223\files_RCH-WA21093-PRIM1\Files\Analysis Ready Datasets\SAS_analysis";
libname r93p2 "E:\SDT-RP-11223\files_RCH-WA21093-PRIM2\Files\Analysis Ready Datasets\SAS_analysis";
libname r93ud "E:\SDT-RP-11223\files_RCH-WA21493-UPDT\Files\Analysis Ready Datasets\SAS_analysis";
libname s223  "E:\SDT-RP-11223\files_SAN-CAMMS223\Files\Datasets_CAMMS223";
libname s323  "E:\SDT-RP-11223\files_SAN-CAMMS323\Files\Datasets_CAMMS323";
libname s507  "E:\SDT-RP-11223\files_SAN-CAMMS32400507\Files\Datasets_CAMMS32400507";

%macro basicdesc(study,sheet,ds,vars1,vars2);
/*  study  - Study-Lib-Name
    sheet  - domain (only for file name)
    ds     - Datasetname
    vars1  - Variablenliste mit nummerischen Variablen (PROC MEANS)
    vars2  - Variablenliste mit Kategorie-/Textvariablen (PROC FREQ)
    Ausgabe als RTF-File
*/
    ods rtf file = "D:\Users\mse0jxh8\Documents\My SAS Files\9.4\basic\&study._&sheet._&ds..rtf" style=styles.journal bodytitle startpage = YES;
    options nodate;
    %let var_list = &vars1.
                    ;
	%let b=1;
	%do %while(%scan(&var_list.,&b.) ne );
	    %let y = %scan(&var_list.,&b.);

		title1 height = 16pt  bold "Variable &y.";

		data ohneformat;
		    set &study..&ds.;
			format _all_ ;
		run;

        proc means data = ohneformat n min q1 mean median q3 max;
            var &y.;
        run;

	    %let b = %eval(&b. + 1);
	%end;

    %let var_list = &vars2.
                    ;
	%let b=1;
	%do %while(%scan(&var_list.,&b.) ne );
	    %let y = %scan(&var_list.,&b.);

		title1 height = 16pt  bold "Variable &y.";

        proc freq data = ohneformat;
            tables &y. /MISSING;
        run;

	    %let b = %eval(&b. + 1);
	%end;

	ods rtf close;
%mend basicdesc;

%basicdesc(N2301,QoL,A_QOL,ACT4C ANX1C DAYNI_1N MOB1C PAI4C QOLTYP1C SCOREBL SFC1C THMSCO1N VASBL,
                     ACT4C ANX1C DAYNI_1N MOB1C PAI4C QOLTYP1C SCOREBL SFC1C STYSID1A THMSCO1N VASBL);
%basicdesc(N2301,OPH-VFT,A_OPH,DAYNI_1N OPHBASEN VADEC1N VAMEAS1C VAOTT1C VATESTC,
                     DAYNI_1N OPHBASEN STYSID1A VADEC1N VAMEAS1C VAOTT1C VASNE1A VATESTC);
%basicdesc(N2301,MSFC,A_MSFC,DAYNI_1N FLAGMSF HPTBL PASATBL  T25FWBL,
                     DAYNI_1N FLAGMSF HPTBL PASATBL STYSID1A T25FWBL);


%basicdesc(N2301,Others,A_MSHIS,DGN_DMS DURMS1 MSLSYM1C MSUP1Y NUMRLP1N NUMRLP2N PEY1C PRVTMT3C 
                                PRVTMT8C RCE5C RELTIME SEV2C SEX1C SYMPRT1C,
                                DGN_DMS DURMS1 MSLSYM1C MSUP1Y NDGCOD1A NUMRLP1N NUMRLP2N PEY1C PRVTMT3C 
                                PRVTMT8C RCE5C RELTIME SEV2C SEX1C SID1A SYMPRT1C NDGCOD1A PT_TXT);


%basicdesc(N2301,AE,A_AEINF,ACNTAK3N AEVSER1C ANALFLAG DAYNI_1N INFTYP2C MHLGCODE SBJ1N VIS1N,
                     ACNTAK3N AEVSER1C ANALFLAG DAYNI_1N INFTYP2C MHLGCODE NMLCOD1A ORGNAM1A SBJ1N SOC_ABBR STY1A VIS1N VISNAM1A);
%basicdesc(N2301,CoMed,A_CMDATC,CONCTIM VIS1N,
                     ATCCODE ATC_TXT CONCTIM SID1A STY1A VIS1N);
%basicdesc(N2301,CoDis,A_CND,ACTPRB1C MHLGCODE,
                     ACTPRB1C MHLGCODE PSTBL_1C SID1A SOC_ABBR STY1A);
%basicdesc(N2301,Demo,A_DMG,BMI_1N ETH2C HGT_1N RCE5C SEX1C VIS_1N VIS_1O WGTBL_1N,
                     BMI_1N ETH2C HGT_1N RCE5C SEX1C SID1A STY1A VIS_1N VIS_1O WGTBL_1N age_1n_cat agedrv1n_cat);
%basicdesc(N2301,IDENT,A_IDENT,BMI_1N CTR1N DCNRSN1C HGT_1N INSTUDY1 INSTUDY4 ITT PK PP PRIORMS RCE5C RDP SAFETY SEX1C TRTN,
                     ASR BMI_1N COU1A CTR1N DCNRSN1C HGT_1N INSTUDY1 INSTUDY4 ITT PK PP PRIORMS RCE5C RDP SAFETY SEX1C SID1A
                     STY1A TGP1A TGPDSC1A TRTN TRTREG1C age_1n_cat);
%basicdesc(N2301,EDSS,A_EDSS,CENPG3N1 CENPG3N2 CENPG6N1 CODPG3N1 CODPG6N1 CONFEXD3 CONFEXD6 CUT12M2O CUTOFF2O C_CDPG3N
                     C_CDPG6N EXDNCHG EXDSC1N EXDSC1NB FCDPG3A2 FCDPG3N1 FCDPG3N2 FCDPG3NA FCDPG6N1 FSEPRF1C FSEPRF1_
                     FUNBOW1N FUNBRA1N FUNCER1N FUNMEN1N FUNPYR1N FUNSEN1N FUNVSU1N F_CDPG3A F_CDPG3N F_CDPG6N KMCENSOR
                     ONSTEDS1 ONSTEDS2 ONSTEXD3 ONSTEXD6 RPEVIS1N TCDPG3N1 TCDPG3N2 TCDPG6N1 TIMPNT1N
                     TIMPNT2N TIMPNT3N T_CDPG3N T_CDPG6N VIS1N VIS1NITT,
                     CENPG3N1 CENPG3N2 CENPG6N1 CODPG3N1 CODPG6N1 CONFEXD3 CONFEXD6 CUT12M2O CUTOFF2O C_CDPG3N
                     C_CDPG6N EXDNCHG EXDSC1N EXDSC1NB FCDPG3A2 FCDPG3N1 FCDPG3N2 FCDPG3NA FCDPG6N1 FSEPRF1C FSEPRF1_
                     FUNBOW1N FUNBRA1N FUNCER1N FUNMEN1N FUNPYR1N FUNSEN1N FUNVSU1N F_CDPG3A F_CDPG3N F_CDPG6N KMCENSOR
                     ONSTEDS1 ONSTEDS2 ONSTEXD3 ONSTEXD6 PERIOD RPEVIS1N SID1A TCDPG3N1 TCDPG3N2 TCDPG6N1 TIMPNT1N
                     TIMPNT2N TIMPNT3N T_CDPG3N T_CDPG6N VIS1N VIS1NITT);
%basicdesc(N2301,MRI,A_MRI,FLAGMRI FREMRI FRET1ITT FRET2ITT FRET2PP FSTERPHA MET3C NUMGD1BA NUMRLP2N RSLVAL1N SCANITT
                     SCANPP TVOLGDBA TVOLT2BA TVOT2PCH TVT1HIBA VIS1N,
                     FLAGMRI FREMRI FRET1ITT FRET2ITT FRET2PP FSTERPHA MET3C NUMGD1BA NUMRLP2N SCANITT
                     SCANPP SID1A VIS1N);
%basicdesc(N2301,Rlps,A_REL,ARLNUM1 ARLNUM2 ARLNUM1C ARLNUM1I ARLNUMA ARLNUMAC ARLNUMB ARLNUMBC ARR1 ARR1I ARRA
                     ARRB ARR_ALL ARR_ALLI CUT12M1O CUTOFF1O DAY_1N INSTUDY1 INSTUDY2 INSTUDYA INSTUDYB INSTUDYI
                     RELCEN1 RELCEN2 RELNUM1 RELNUM2 RELNUM1C RELNUM1I RELNUMA RELNUMAC RELNUMB RELNUMBC,
                     ARLNUM1 ARLNUM2 ARLNUM1C ARLNUM1I ARLNUMA ARLNUMAC ARLNUMB ARLNUMBC ARR1 ARR1I ARRA
                     ARRB ARR_ALL ARR_ALLI CUT12M1O CUTOFF1O DAY_1N INSTUDY1 INSTUDY2 INSTUDYA INSTUDYB INSTUDYI
                     RELCEN1 RELCEN2 RELNUM1 RELNUM2 RELNUM1C RELNUM1I RELNUMA RELNUMAC RELNUMB RELNUMBC SID1A);
%basicdesc(N2301,Others,A_MSHIS,DGN_DMS DURMS1 MSLSYM1C MSUP1Y NUMRLP1N NUMRLP2N PEY1C PRVTMT3C PRVTMT8C
                     RCE5C RELTIME SEV2C SEX1C SYMPRT1C,
                     DGN_DMS DURMS1 MSLSYM1C MSUP1Y NDGCOD1A NUMRLP1N NUMRLP2N PEY1C PRVTMT3C PRVTMT8C
                     RCE5C RELTIME SEV2C SEX1C SID1A SYMPRT1C);
%basicdesc(N2301,Ophta,A_OPHHIS,NEUBILAT NEURECUR NEURITIS NEUSNEPI NEUUNLAT,
                     NEUBILAT NEURECUR NEURITIS NEUSNEPI NEUUNLAT SID1A);
%basicdesc(N2301,Lab,A_LRS,BLVAL_1N CVLLN_1C CVRSL_1N CVULN_1C CV_BL_1N LABRSL1N VIS1N,
                     BLFLG_1C CVLLN_1C CVULN_1C CVUNT_1C LABCAT1C MTYPE_1C
                     PARNAM1C PARUNT1C PREUNT1C PSTBL_1C SID1A STY1A VIS1N);



