function [rho_H2O, mu_H2O, k_H2O, C_H2O] = water(T_H2O)
%CoolProp water properties at 5 bara as function of T

properties = [...
290	998.9885861	0.001083794	0.592539544	4185.285586
300	996.7358498	0.000853707	0.60972339	4179.510429
310	993.5591948	0.000693365	0.624483374	4178.233451
320	989.6009056	0.000576799	0.637204771	4179.599423
330	984.960843	0.00048924	0.648119352	4182.758744
340	979.7114128	0.000421735	0.657377739	4187.420596
350	973.9061774	0.000368576	0.665087467	4193.597382
360	967.5850682	0.000325964	0.671333053	4201.455891
370	960.7776521	0.000291283	0.676186364	4211.234329
380	953.5052229	0.000262682	0.679711727	4223.200463
390	945.7821457	0.000238817	0.681968212	4237.634951
400	937.6167009	0.000218691	0.68301042	4254.829758
410	929.0115701	0.000201554	0.682888551	4275.095605
420	919.9640437	0.00018683	0.681648122	4298.77518
430	2.629677412	1.42369E-05	0.030979346	2346.867012
440	2.557901804	1.46604E-05	0.031804977	2267.13806
450	2.491282417	1.5083E-05	0.032645223	2215.892459
460	2.429002993	1.55051E-05	0.033500427	2178.712819
470	2.370491023	1.59268E-05	0.034370882	2150.424842
480	2.315304325	1.63481E-05	0.035256736	2128.541336
490	2.263086111	1.67692E-05	0.036158009	2111.550712
500	2.213541394	1.71902E-05	0.03707461	2098.404618
510	2.16642211	1.7611E-05	0.038006369	2088.328909
520	2.12151682	1.80316E-05	0.03895305	2080.73424
530	2.078643227	1.84521E-05	0.039914376	2075.164107
540	2.037642594	1.88724E-05	0.040890036	2071.260639
550	1.998375499	1.92926E-05	0.041879695	2068.740555
560	1.960718537	1.97125E-05	0.042883009	2067.377706
570	1.924561746	2.01323E-05	0.043899621	2066.990134
580	1.889806555	2.05519E-05	0.044929175	2067.43037
590	1.856364126	2.09712E-05	0.045971314	2068.578064
600	1.824154018	2.13902E-05	0.047025684	2070.334358
610	1.79310308	2.18089E-05	0.048091938	2072.61753
620	1.763144535	2.22273E-05	0.049169735	2075.359619
630	1.734217212	2.26453E-05	0.050258741	2078.503779
640	1.7062649	2.30628E-05	0.051358634	2082.002183
650	1.679235798	2.348E-05	0.052469097	2085.814377
660	1.653082045	2.38967E-05	0.053589826	2089.905952
670	1.627759311	2.43129E-05	0.054720525	2094.247482
680	1.603226446	2.47286E-05	0.055860908	2098.81367
690	1.579445176	2.51437E-05	0.057010698	2103.582649
700	1.556379831	2.55582E-05	0.058169627	2108.535415
710	1.533997111	2.59722E-05	0.059337438	2113.655365
720	1.512265876	2.63855E-05	0.06051388	2118.927911
730	1.491156961	2.67982E-05	0.061698714	2124.340166
740	1.470643011	2.72102E-05	0.062891706	2129.880681
750	1.450698332	2.76215E-05	0.064092631	2135.53923
760	1.431298762	2.80321E-05	0.065301274	2141.306625
770	1.412421549	2.8442E-05	0.066517424	2147.174565
780	1.394045243	2.88511E-05	0.067740878	2153.135511
790	1.376149602	2.92595E-05	0.068971442	2159.182573
800	1.358715503	2.9667E-05	0.070208925	2165.309429
810	1.341724862	3.00738E-05	0.071453146	2171.510239
820	1.325160557	3.04797E-05	0.072703926	2177.779587
830	1.30900637	3.08848E-05	0.073961095	2184.112427
840	1.29324692	3.12891E-05	0.075224487	2190.504031
850	1.277867606	3.16925E-05	0.076493942	2196.949957
860	1.262854563	3.2095E-05	0.077769303	2203.446011
870	1.248194609	3.24966E-05	0.07905042	2209.988223
880	1.233875203	3.28973E-05	0.080337147	2216.572817
890	1.219884407	3.32971E-05	0.081629342	2223.196199
900	1.206210847	3.3696E-05	0.082926868	2229.854933
910	1.192843681	3.4094E-05	0.08422959	2236.545728
920	1.179772564	3.4491E-05	0.08553738	2243.265428
930	1.166987624	3.48871E-05	0.086850111	2250.010999
940	1.154479431	3.52823E-05	0.088167661	2256.779522
950	1.142238971	3.56764E-05	0.089489911	2263.56818
960	1.130257628	3.60697E-05	0.090816745	2270.374261
970	1.118527156	3.64619E-05	0.092148046	2277.195142
980	1.10703966	3.68532E-05	0.093483716	2284.028291
990	1.095787582	3.72435E-05	0.094823645	2290.871262
1000	1.084763677	3.76328E-05	0.096167726	2297.721689
1010	1.073960997	3.80211E-05	0.097515858	2304.577283
1020	1.06337288	3.84085E-05	0.098867942	2311.435833
1030	1.052992933	3.87948E-05	0.100223883	2318.295198
1040	1.042815016	3.91802E-05	0.101583587	2325.153312
1050	1.032833232	3.95645E-05	0.102946962	2332.008175
1060	1.023041914	3.99479E-05	0.104313921	2338.857856
1070	1.013435614	4.03302E-05	0.105684375	2345.700488
1080	1.004009094	4.07116E-05	0.10705824	2352.534271
1090	0.99475731	4.10919E-05	0.108435433	2359.357467
1100	0.985675412	4.14712E-05	0.109815874	2366.168399
1110	0.976758725	4.18496E-05	0.111199483	2372.965453
1120	0.96800275	4.22269E-05	0.112586184	2379.747071
1130	0.959403148	4.26032E-05	0.1139759	2386.511757
1140	0.950955739	4.29785E-05	0.11536856	2393.25807
1150	0.94265649	4.33528E-05	0.116764089	2399.984626
1160	0.934501512	4.37261E-05	0.118162419	2406.690096
1170	0.92648705	4.40984E-05	0.11956348	2413.373206
1180	0.918609483	4.44697E-05	0.120967205	2420.032734
1190	0.910865311	4.48399E-05	0.122373528	2426.66751
1200	0.903251154	4.52092E-05	0.123782384	2433.276414
];

tempinterp = properties(:,1);
interp = interp1q(tempinterp, properties, T_H2O);

rho_H2O = interp(2); %density, kg/m^3
mu_H2O = interp(3); %viscocity, Pa*s
k_H2O = interp(4); %conductive heat trans coef (W/mK)
C_H2O = interp(5); %specific heat capacity (J/kgK)
