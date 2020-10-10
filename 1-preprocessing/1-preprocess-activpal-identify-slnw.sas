/*=============================================*/
/*=================REFERENCE===================*/
/*=============================================*/
/*Source: 
Winkler, E. A., Bodicoat, D. H., Healy, G. N., Bakrania, K., Yates, T., 
Owen, N., Dunstan, D. W., & Edwardson, C. L. (2016).
Identifying adults' valid waking wear time by automated estimation in activPAL data collected with a 24 h wear protocol. 
Physiological measurement, 37(10), 1653–1668. https://doi.org/10.1088/0967-3334/37/10/1653*/

*Define file path;
%let FILEPATH = "/folders/myfolders/dataset/activpal_pid_61_72_fu.csv"; 
%let PROCESSED_FILEPATH = "/folders/myfolders/dataset/sas-processed/activpal_pid_61_72_fu_processed.csv";

proc import datafile = &FILEPATH dbms=csv out=activpal replace;
	guessingrows=max;
run;

proc contents data=activpal;
run;

/*====================================*/
/*=========FORMAT INPUT TYPE==========*/
/*====================================*/

data activpal_format;
set WORK.activpal; 

format time best32.; 
format interval_s best32. ; 
format activity_code 1.0 ; 
format cumulativestepcount best32. ; 
format activity_score_met_h best32.;
format datacount_samples best12. ;
format patient_id best8.;
keep time datacount_samples interval_s
 	activity_code
 	cumulativestepcount activity_score_met_h 
 	patient_id visit_info;
run;

/**/
data activpal_format; 
set activpal_format; 
if Time in (0,.) 
	then delete; *remove rows of 0s from end of some activpal event files;
boutid = _N_; 
sasdatetimenumber = (Time- 21916)*86400; *convert exceldate to sasdatetime;
cumulativestepcount_i = cumulativestepcount ; 
steps_n = cumulativestepcount - lag(cumulativestepcount_i); 
if boutid = 1 then steps_n = cumulativestepcount ; 
drop cumulativestepcount_i ; *non cumulative steps;
run;

*recalculate intervals as time between each bout -> more precise than the rounded data;
proc sort data=activpal_format; 
by patient_id descending boutid; 
run;

data activpal_format; 
set activpal_format;
boutidR = _N_;
sasdatetimenumber_endbout = lag(sasdatetimenumber); 
if sasdatetimenumber_endbout = . & sasdatetimenumber ne . 
	then sasdatetimenumber_endbout = sasdatetimenumber + interval_s; 
interval_s_derived = sasdatetimenumber_endbout - sasdatetimenumber;
basicdate = datepart(sasdatetimenumber); 
noondate = datepart(sasdatetimenumber + 43200); *dates at beginning of bout for day and for noon to noon window;
format sasdatetimenumber sasdatetimenumber_endbout datetime. basicdate noondate ddmmyy10.;
run;

/*  */
*PART 2: The algorithm steps;
proc sort data = activpal_format; 
by patient_id basicdate boutid; 
run;

/*
identifying intial sleep: 
longest bout >= &longest at least (e.g. 2h) that occurred in each 24-hr window (defined noon-noon) 
and any bouts >= &vlongbout (e.g. 5h)
*/
proc means data = activpal_format noprint; 
by patient_id noondate; 
output out = noondate max(interval_s_derived)=maxboutdur_nd; 
run;

%let longestatleast = 2;
%let vlongbout = 5;
%let longbout = 2;
%let mediumbout = 0.5;

data activpal_format; 
merge activpal_format noondate; 
by patient_id noondate;
if interval_s_derived = maxboutdur_nd & 
	(interval_s_derived / (60*60)) >= &longestatleast 
	then longestbout = 1;
sleepbout = 0; 
if longestbout = 1 
	then sleepbout = 1; 
if (interval_s_derived/(60*60)) ge &vlongbout 
	then sleepbout = 1;
	
*initial sleep yes/no;
*calculating the various times in each activity for summation (to identify probable nonwear days) and for searching
around sleep;
sitting_s = 0; 
standing_s = 0; 
stepping_s = 0; 
active_s = 0; 
nonstep_s = 0; 
prnonstep2h_s = 0; 
prnonstep30m_s = 0;
if activity_code = 0 
	then sitting_s = interval_s_derived; 
if activity_code = 1 
	then standing_s = interval_s_derived; 
if activity_code = 2 
	then stepping_s = interval_s_derived;
if activity_code in (1,2) 
	then active_s = interval_s_derived; 
if activity_code in (0,1) 
	then nonstep_s = interval_s_derived;
if activity_code in (0,1) & 
	(interval_s_derived/(60*60)) ge &longbout 
	then prnonstep2h_s = interval_s_derived;
if activity_code in (0,1) & 
	(interval_s_derived/(60*60)) ge &mediumbout 
	then prnonstep30m_s = interval_s_derived;
drop _type_ _freq_;
run;

/*  
*numbering the periods between one inital sleep bout and the next;
*/
data activpal_format; 
set activpal_format; 
by patient_id;
if first.patient_id 
	then slnum = 0; 
if sleepbout = 1 
	then slnum + 1; *periods forwards;
sleepboutall = 0; 
if sleepbout = 1 
	then sleepboutall = sleepbout; *sleepboutall refers to all sleep including newly found;
slnumR = -1*(slnum); 
slnumR_i = lag(slnumR); 
if sleepbout = 1 
	then slnumR =slnumR_i ; 
	drop slnumR_i; *periods backwards;
run;

*ALGORITHM STEP 2 -> SEARCHING AROUND INTIAL SLEEP, AFTER (FORWARDS) AND BEFORE (BACKWARDS);
*this version starts at the first bout, continues whether more sleep is not found or not, then later the windows are
limited to the relevant periods
- only more sleep following each previously identified sleep bout until more sleep is not found (once) ;
*FORWARDS;

%let steptolerance = 20;
%let sleepwindow = 15;
data psp; 
	set activpal_format; 
	by patient_id boutid ;
	if first.patient_id then ps_num=0; *post sleep window period number;
	if first.patient_id or moresleep or notmoresleep 
		then do;
		start=0; *indicator for starting to count the sleep period;
		moresleep = 0; *indicator for stopping the window - and it was sleep;
		notmoresleep = 0; *indicator for stopping the window - and it was not more sleep;
		cnt_steps=0; *number of steps in window;
		dur_step_s=0; *duration of stepping in window;
		dur_nonstep_s=0; *duration of sitting or standing in window;
		dur_prnonstep2h_s = 0; *duration of 2+hr bouts of sitting or standing in window;
		dur_prnonstep30m_s = 0;*duration of 30+min bouts of sitting or standing in window;
		n_sleepbouts = 0; *number of previously identifed sleep bouts found in window;
		strt_ps_num = boutid; *bout number at start of window;
		end_ps_num = boutid; *bout number at end of window;
		strt_ps_sb = sasdatetimenumber; *start of window, bout start time;
		end_ps_sb = sasdatetimenumber; *start of window, bout end time;
		strt_ps_eb = sasdatetimenumber_endbout; *end of window, bout start time;
		end_ps_eb = sasdatetimenumber_endbout; *end of window, bout end time;
		s_window = 0 ; *window duration [after sleep/last window to beginning of bout at end of this window];
	end;
	retain ps_num start moresleep notmoresleep strt_ps_eb strt_ps_sb end_ps_eb end_ps_sb s_window cnt_steps
		dur_step_s dur_nonstep_s dur_prnonstep2h_s dur_prnonstep30m_s n_sleepbouts strt_ps_num end_ps_num boutid ;

	*The beginning of the window (initially, or after a prior iteration);
	if start = 0 
		then do; 
		strt_ps_eb=sasdatetimenumber_endbout; 
		strt_ps_sb=sasdatetimenumber; 
		strt_ps_num = boutid;
		start=1; 
		end;
	
	*keep track of times and bout numbers and count the steps, sleep bouts durations etc as the window continues until it
	ends;
	end_ps_sb=sasdatetimenumber; 
	end_ps_eb=sasdatetimenumber_endbout; 
	end_ps_num = boutid;
	cnt_steps + steps_n; 
	dur_step_s + stepping_s; 
	dur_sit_s + sitting_s; 
	dur_nonstep_s + nonstep_s;
	dur_prnonstep2h_s + prnonstep2h_s ; 
	dur_prnonstep30m_s + prnonstep30m_s ;
	if sleepboutall = 1 
		then n_sleepbouts = n_sleepbouts + 1;
	s_window = end_ps_sb - strt_ps_eb;

	*STOP the window and begin another one if conditions for more sleep are met or the window ends without meeting the
	conditions;
	*STOPPED: found another sleep bout or a bout >= &longbout (eg 2h) within 15 mins (&sleepwindow). -> moresleep.;
	*STOPPED: found a bout >= &mediumbout (e.g. >=0.5h) within (&sleepwindow). -> moresleep if found < 20 steps first
	otherwise notmoresleep;
	*STOPPED: scanned all window (&sleepwindow) eg 15 min or the last bout found. -> moresleep if found 0 steps;
	if
		( (s_window/60) <= &sleepwindow & (n_sleepbouts > 0 | dur_prnonstep2h_s > 0)) |
		( (s_window/60) <= &sleepwindow & (dur_prnonstep30m_s > 0 | n_sleepbouts > 0))|
		( (s_window/60) >= &sleepwindow) | last.patient_id
		then do;
	* found condition for more sleep;
	if ((s_window/60) <= &sleepwindow & (n_sleepbouts > 0 | dur_prnonstep2h_s > 0)) |
		((s_window/60) <= &sleepwindow & (dur_prnonstep30m_s > 0 | n_sleepbouts > 0) & cnt_steps < &steptolerance ) |
		cnt_steps = 0
		then moresleep=1;
	*met no conditions for finding more sleep;
	else 
		moresleep=0; 
	notmoresleep = 1- moresleep;
	ps_num = ps_num+1; *number the potential sleep period;
	
	*output one record for each potential sleep period;
	keep patient_id ps_num strt_ps_eb strt_ps_sb end_ps_eb end_ps_sb s_window moresleep notmoresleep cnt_steps dur_step_s
	dur_nonstep_s dur_prnonstep2h_s dur_prnonstep30m_s n_sleepbouts strt_ps_num end_ps_num ;
	output;
	end;
format strt_ps_eb strt_ps_sb end_ps_eb end_ps_sb datetime.;
run;

/*  
*decide whether each window is extra sleep or not -> only extra sleep if found more sleep but not if searching forwards
before found any sleep at all;
*/
data psp;
	set psp; 
	by patient_id; 
	if first.patient_id 
		then periodnum = 0; 
	if lag(n_sleepbouts) > 0 
		then periodnum+1; 
	if notmoresleep = 1
		then stopnum1 = ps_num; 
run;

proc means data = psp noprint; 
by patient_id periodnum; 
where notmoresleep= 1; 
output out=stops min(stopnum1)=firststop1; 
run;

data psp; 
merge psp stops; 
by patient_id periodnum; 
extrasleep1 = moresleep; 
if periodnum = 0 then extrasleep1 = 0; 
if ps_num ge firststop1 & firststop1 ne . 
	then extrasleep1 = 0; 
run;

data psp1; 
	set psp; 
	if periodnum ne 0 & extrasleep1 = 1; 
run;

*list out all the new sleep bouts found sarching forwards and mark those bouts as yes to sleepboutall;
data psplist(keep=patient_id boutid sleepboutF); 
	set psp1; 
	by patient_id ps_num;
	do i=strt_ps_num to end_ps_num by 1; 
	boutid = i; 
	sleepboutF=1; 
	output; 
	end;
run;

data activpal_format; 
	merge activpal_format(in = in_events) psplist; 
	by patient_id boutid; 
	if in_events; 
	if (sleepboutF = 1 | sleepbout = 1)
		then sleepboutall = 1; 
run;

/*
* LARGELY THE SAME, BACKWARDS;
*/

proc sort data = activpal_format; 
by patient_id boutidR; 
run;

data pspR; 
	set activpal_format; 
	by patient_id boutidR ; *DIFFERENT BACKWARDS;
	if first.patient_id 
		then ps_num=0;
	if first.patient_id or moresleep or notmoresleep 
		then do;
			start=0; 
			moresleep = 0; 
			notmoresleep = 0; 
			cnt_steps=0; 
			dur_step_s=0; 
			dur_nonstep_s=0;
			dur_prnonstep2h_s = 0; 
			dur_prnonstep30m_s = 0; 
			n_sleepbouts = 0;
			strt_ps_num = boutidR; *DIFFERENT BACKWARDS;
			end_ps_num = boutidR; *DIFFERENT BACKWARDS;
			strt_ps_sb = sasdatetimenumber; 
			end_ps_sb = sasdatetimenumber; 
			strt_ps_eb = sasdatetimenumber_endbout;
			end_ps_eb = sasdatetimenumber_endbout; 
			s_window = 0 ;
		end;
	retain ps_num start moresleep notmoresleep strt_ps_eb strt_ps_sb end_ps_eb end_ps_sb s_window cnt_steps dur_step_s
		dur_nonstep_s dur_prnonstep2h_s dur_prnonstep30m_s n_sleepbouts strt_ps_num end_ps_num boutid ;
	if start = 0 
		then do; 
		strt_ps_eb=sasdatetimenumber_endbout; 
		strt_ps_sb=sasdatetimenumber; 
		strt_ps_num = boutidR;
		start=1; 
	end; *DIFFERENT BACKWARDS;
	end_ps_sb=sasdatetimenumber; 
	end_ps_eb=sasdatetimenumber_endbout; 
	end_ps_num = boutidR; *DIFFERENT BACKWARDS;
	cnt_steps + steps_n; 
	dur_step_s + stepping_s; 
	dur_sit_s + sitting_s; 
	dur_nonstep_s + nonstep_s;
	dur_prnonstep2h_s + prnonstep2h_s ; 
	dur_prnonstep30m_s + prnonstep30m_s ;
	if sleepboutall = 1 
		then n_sleepbouts = n_sleepbouts + 1;
		s_window = (strt_ps_sb - end_ps_eb); *DIFFERENT BACKWARDS: time since beginning of period, with window running backwards;
	if
		((s_window/60) <= &sleepwindow & (n_sleepbouts > 0 | dur_prnonstep2h_s > 0)) |
		((s_window/60) <= &sleepwindow & (dur_prnonstep30m_s > 0 | n_sleepbouts > 0))|
		((s_window/60) >= &sleepwindow) |
		last.patient_id
	then do;
		if ((s_window/60) <= &sleepwindow & (n_sleepbouts > 0 | dur_prnonstep2h_s > 0)) |
			((s_window/60) <= &sleepwindow & (dur_prnonstep30m_s > 0 | n_sleepbouts > 0) & cnt_steps < &steptolerance) |
			cnt_steps = 0
		then moresleep=1;
		else moresleep=0; 
	notmoresleep = 1- moresleep;
	ps_num=ps_num+1;
	keep patient_id ps_num strt_ps_eb strt_ps_sb end_ps_eb end_ps_sb s_window moresleep notmoresleep cnt_steps
	dur_step_s dur_nonstep_s dur_prnonstep2h_s dur_prnonstep30m_s n_sleepbouts strt_ps_num end_ps_num ;
	output;
	end;
format strt_ps_eb strt_ps_sb end_ps_eb end_ps_sb datetime.;
run;

*decide whether each window is extra sleep or not -> only extra sleep if found following (backwards) from initial sleep
or moresleep without first encountering a failure to find moresleep;
data pspR; 
	set pspR; 
	by patient_id;
	if first.patient_id 
		then periodnumR = 0; 
	if lag(n_sleepbouts) > 0 then periodnumR+1; 
	if notmoresleep = 1 then stopnum1 = ps_num; 
run;

proc means data = pspR noprint; 
by patient_id periodnumR; 
output out=stopsR min(stopnum1)=firststop1 ; 
run;

data pspR; 
merge pspR stopsR; 
by patient_id periodnumR; 
extrasleep1 = moresleep; 
if periodnumR = 0 then extrasleep1 = 0; 
if ps_num ge firststop1 & firststop1 ne . 
	then extrasleep1 = 0; 
run;

data pspR1; 
	set pspR; 
	if periodnumR ne 0 & extrasleep1 = 1; 
run;

data psplistR(keep=patient_id boutidR sleepboutB); 
	set pspR1; 
	by patient_id ps_num;
	do i=strt_ps_num to end_ps_num by 1; 
		boutidR=i; 
		sleepboutB=1; 
		output; 
	end;
run;

proc sort data = psplistR; 
by patient_id boutidR; 
run;

proc sort data = activpal_format; 
by patient_id boutidR; 
run;

data activpal_format; 
merge activpal_format(in = in_events) psplistR; 
by patient_id boutidR; 
if in_events;
if (sleepboutF = 1 | sleepboutB = 1 | sleepbout = 1) 
	then sleepboutall = 1;
run;

proc sort data = activpal_format; 
by patient_id basicdate boutid; 
run;

*daily summary to determine valid days -> durations in mins, and values of 0 set to 0 not missing, and percentages of
waking wear;
data activpal_format; 
	set activpal_format;
	*time in each activity during waking wear;
	AW_sitting_t = 0; 
	AW_standing_t = 0; 
	AW_stepping_t = 0; 
	AW_active_t = 0; 
	AW_all_t = 0; 
	AW_steps_n = 0;
	if sleepboutall ne 1
		then do; 
		AW_sitting_t = sitting_s; 
		AW_standing_t = standing_s; *AW_standing_t = sitting_s; 
		AW_stepping_t = stepping_s;
		AW_active_t = active_s; 
		AW_all_t = interval_s_derived; 
		AW_steps_n = steps_n; 
		end;
	*time in each activity during sleep or nonwear;
	SL_sitting_t = 0 ; 
	SL_standing_t = 0; 
	SL_stepping_t = 0 ; 
	SL_active_t = 0 ; 
	SL_all_t = 0; 
	SL_steps_n = 0;
	if sleepboutall = 1 
		then do; 
		SL_sitting_t = sitting_s; 
		SL_standing_t = standing_s; *SL_standing_t = sitting_s; 
		SL_stepping_t = stepping_s;
		SL_active_t = active_s; 
		SL_all_t = interval_s_derived; 
		SL_steps_n = steps_n; 
		end;
run;

proc means data = activpal_format noprint; 
by patient_id basicdate; 
output out = summary
	sum(AW_sitting_t)=AW_sitting_t_TOT 
	sum(AW_standing_t)=AW_standing_t_TOT 
	sum(AW_stepping_t)=AW_stepping_t_TOT
	sum(AW_active_t)=AW_active_t_TOT 
	sum(AW_all_t)=AW_all_t_TOT 
	sum(AW_steps_n)=AW_steps_n
	sum(SL_sitting_t)=SL_sitting_t_TOT 
	sum(SL_standing_t)=SL_standing_t_TOT 
	sum(SL_stepping_t)=SL_stepping_t_TOT
	sum(SL_active_t)=SL_active_t_TOT 
	sum(SL_all_t)=SL_all_t_TOT 
	sum(SL_steps_n)=SL_steps_n;
run;

%let minhrsvalid = 10;
%let minstepsday = 500;
%let maxdaypct = 95;

data summary; 
set summary; 
drop _freq_ _type_;
array durs SL_sitting_t_TOT SL_standing_t_TOT SL_stepping_t_TOT SL_active_t_TOT SL_all_t_TOT AW_sitting_t_TOT
	AW_standing_t_TOT AW_stepping_t_TOT AW_active_t_TOT AW_all_t_TOT ;
do over durs; 
	if durs = . then durs = 0; 
	durs = durs / 60; 
end;
if AW_all_t_TOT > 0 
	then do;
	AW_sitting_t_TOTPCT = 100*(AW_sitting_t_TOT / AW_all_t_TOT); 
	AW_standing_t_TOTPCT = 100*(AW_standing_t_TOT / AW_all_t_TOT); 
	AW_stepping_t_TOTPCT =100*(AW_stepping_t_TOT / AW_all_t_TOT);
	AW_active_t_TOTPCT =100*(AW_active_t_TOT / AW_all_t_TOT);
end;
array oth SL_steps_n AW_steps_n AW_sitting_t_TOTPCT AW_standing_t_TOTPCT AW_stepping_t_TOTPCT AW_active_t_TOTPCT;
do over oth; 
	if oth = . then oth = 0; 
end;
*indicate whether day was invalid wear day;
validday = 1;
if AW_all_t_TOT < 60*&minhrsvalid 
	then validday = 0; *not valid if <= &minhrsvalid waking wear (e.g. 10h);
if (AW_standing_t_TOTPCT ge &maxdaypct | AW_sitting_t_TOTPCT ge &maxdaypct | AW_stepping_t_TOTPCT ge &maxdaypct) 
	then validday = 0; *not valid if >= &maxdaypct (e.g. 95%) in one activity;
if AW_steps_n < &minstepsday 
	then validday = 0; *not valid if < &minstepsday (e.g. 500);
format basicdate ddmmyy10.;
run;

data activpal_format; 
	merge activpal_format summary; 
	by patient_id basicdate; 
run;

/* %let PROCESSED_FILEPATH = "/folders/myfolders/dataset/activpal_pid_1_2_processed.csv"; */
proc export data = activpal_format
	outfile = &PROCESSED_FILEPATH dbms=csv replace;
	putnames = yes;
run;

/*  */
/*  */
/* *PART 3: Outputting the data to csv for later use or quality control steps; */
/*  */
/* %let csvname = /folders/myfolders/dataset/activpal_pid_1_2_sas; */
/* %let tagfiles = 1; */
/*  */
/* %if &tagfiles = 1 %then %do; */
/* *output csv files like the events information and new sleep/nonwear and valid day classifications; */
/* 	data taggedevents;  */
/* 		set activpal_format; */
/* 		APDatetimevar = (sasdatetimenumber/86400) + 21916 ;  */
/* 		APDatacount = datacount_samples;  */
/* 		APInterval_s = interval_s; */
/* 		APActivity_code = activity_code;  */
/* 		APCumulativeStepCount = CumulativeStepCount;  */
/* 		APActivityScoreMEThr = activityscore_met_h; */
/* 		validday_est = validday;  */
/* 		sleepbout_est = sleepboutall;  */
/* 		participantID = patient_id; */
/* 		keep APDatetimevar APDatacount APInterval_s APActivity_code APCumulativeStepCount APActivityScoreMEThr validday_est */
/* 		sleepbout_est participantID; */
/* 	run; */
/*  */
/* 	proc export data = taggedevents  */
/* 	outfile = "&csvname..csv" dbms=csv replace; putnames = yes; run; */
/* 	%end; */
/*  */
/* %let outlist = 1; */
/* %let listname = /folders/myfolders/dataset/activpal_pid_1_2_list; */
/* %if &outlist = 1 %then %do; */
/* 	*outputs a csv list of the sleep-nonwear periods [SLNW bouts] for that id; */
/* 	data forlist (keep = patient_id slbnum sasdatetimenumber sasdatetimenumber_endbout sleepboutall);  */
/* 		set activpal_format;  */
/* 		by patient_id boutid; */
/* 		if first.patient_id then do;  */
/* 			slbnum = 1;  */
/* 			end; */
/* 		if sleepboutall ne lag(sleepboutall) then do;  */
/* 			slbnum + 1;  */
/* 		end; */
/* 	run; */
/* 	proc means data = forlist noprint;  */
/* 	by patient_id slbnum;  */
/* 	where sleepboutall = 1;  */
/* 	output out=listsl min(sasdatetimenumber)=slnwperiod_st max(sasdatetimenumber_endbout)=slnwperiod_end;  */
/* 	run; */
/*  */
/* 	data listsl;  */
/* 		set listsl;  */
/* 		drop _type_ _freq_;  */
/* 		format slnwperiod_st slnwperiod_end datetime.;  */
/* 		drop slbnum;  */
/* 	run; */
/* 	proc export data= listsl outfile= "&listname..csv" dbms=csv replace; putnames=yes;  */
/* 	run; */
/* 	%end; */
/* 	 */
/* %let keeplist = 1; */
/* %let keepname = /folders/myfolders/dataset/activpal_pid_1_2_keep_list; */
/* %if &keeplist = 1 %then %do; */
/* 	*outputs a csv list of all the periods of data to keep as valid waking wear [not SLNW bouts, not invalid days] for that */
/* 	id; */
/* 	data forkplist (keep = patient_id valnum sasdatetimenumber sasdatetimenumber_endbout incdata sleepboutall validday);  */
/* 		set activpal_format;  */
/* 		by patient_id boutid; */
/* 		incdata = 1;  */
/* 		if sleepboutall = 1 | validday ne 1  */
/* 			then incdata = 0; */
/* 		if first.patient_id then do;  */
/* 			valnum = 1;  */
/* 		end; */
/* 		if incdata ne lag(incdata) then do;  */
/* 			valnum + 1;  */
/* 		end; */
/* 	run; */
/* 	 */
/* 	proc means data = forkplist noprint;  */
/* 		by patient_id valnum;  */
/* 		where incdata = 1;  */
/* 		output out=listval min(sasdatetimenumber)=validperiod_st max(sasdatetimenumber_endbout)=validperiod_end;  */
/* 	run; */
/* 	 */
/* 	data listval;  */
/* 		set listval;  */
/* 		drop _type_ _freq_;  */
/* 		format validperiod_st validperiod_end datetime.;  */
/* 		drop valnum;  */
/* 	run; */
/* 	 */
/* 	proc export data= listval outfile= "&keepname..csv" dbms=csv replace; putnames=yes;  */
/* 	run; */
/* %end; */