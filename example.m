t = edu.buaa.ACTMapMatching
t.load('D:\map-matching-matlab\cache', 20); % 20 is the gps accuracy, unit one meter
% example trajectory, first column is latitude, second is longitude, third is timestamp (seconds).
trajectory=[
39.862180	116.471353	1478793600
39.862180	116.471353	1478793606
39.862180	116.471353	1478793612
39.862180	116.471353	1478793618
39.862180	116.471353	1478793624
39.862180	116.471353	1478793630
39.862180	116.471353	1478793636
39.862180	116.471353	1478793642
39.862180	116.471353	1478793648
39.862180	116.471353	1478793654
39.862180	116.471353	1478793660
39.862180	116.471353	1478793666
39.862180	116.471353	1478793672
39.862180	116.471353	1478793678
39.862180	116.471353	1478793684
39.862180	116.471353	1478793690
39.862180	116.471353	1478793696
39.862180	116.471353	1478793702
39.862180	116.471353	1478793708
39.862180	116.471353	1478793714
39.862180	116.471353	1478793720
39.862180	116.471353	1478793726
39.862155	116.471312	1478793732
39.862398	116.471278	1478793738
39.862892	116.471308	1478793744
39.862892	116.471308	1478793750
39.863462	116.471302	1478793756
39.864285	116.47125	1478793762
39.864728	116.471197	1478793768
39.864728	116.471197	1478793774
39.865840	116.47124	1478793780
39.865840	116.47124	1478793786
39.866227	116.471272	1478793792
39.866227	116.471272	1478793798
39.866227	116.471272	1478793804
39.866227	116.471272	1478793810
39.866227	116.471272	1478793816
39.866227	116.471272	1478793822
39.866227	116.471272	1478793828
39.866227	116.471272	1478793834
39.866227	116.471272	1478793840
39.866683	116.471233	1478793846
39.867257	116.471232	1478793852
39.867898	116.471197	1478793858
39.868600	116.471212	1478793864
39.869195	116.47123	1478793870
39.869748	116.471265	1478793876
39.869902	116.471512	1478793882
39.869905	116.472032	1478793888
39.869910	116.47254	1478793894
39.869925	116.473153	1478793900
39.869925	116.473153	1478793906
39.869925	116.473153	1478793912
39.869925	116.473153	1478793918
39.869925	116.473153	1478793924
39.869925	116.473153	1478793930
39.869925	116.473153	1478793936
39.869925	116.473153	1478793942
39.869925	116.473153	1478793948
39.869925	116.473153	1478793954
39.869925	116.473153	1478793960
39.869925	116.473153	1478793966
39.869862	116.474017	1478793972
39.869892	116.47471	1478793978
39.869907	116.475665	1478793984
39.869892	116.476687	1478793990
39.869900	116.477612	1478793996
39.869877	116.478773	1478794002
39.869883	116.480057	1478794008
39.869875	116.481368	1478794014
39.869862	116.482427	1478794020
39.869820	116.483377	1478794026
39.869718	116.483988	1478794032
39.869177	116.484415	1478794038
39.868713	116.483857	1478794044
39.868807	116.48328	1478794050
39.869247	116.482982	1478794056
39.870002	116.483012	1478794062
39.871343	116.483233	1478794068
39.872083	116.48327	1478794074
39.872083	116.48327	1478794080
39.873123	116.483387	1478794086
39.874967	116.48359	1478794092
39.876015	116.483678	1478794098
39.877113	116.48371	1478794104
39.878237	116.483737	1478794110
39.878237	116.483737	1478794116
39.880498	116.48374	1478794122
39.881485	116.483742	1478794128
39.882435	116.483763	1478794134
39.883353	116.483762	1478794140
39.884500	116.483765	1478794146
39.885428	116.483772	1478794152
39.886532	116.483778	1478794158
39.887510	116.483793	1478794164
39.888980	116.483797	1478794170
39.890307	116.483787	1478794176
39.891385	116.483798	1478794182
39.892432	116.483797	1478794188
39.893398	116.483798	1478794194
39.894527	116.483808	1478794200
39.896623	116.483795	1478794206
39.897798	116.483785	1478794212
39.899008	116.483792	1478794218
39.899008	116.483792	1478794224
39.900500	116.483793	1478794230
39.901608	116.483813	1478794236
39.902650	116.483807	1478794242
39.904723	116.48382	1478794248
39.905773	116.48381	1478794254
39.906687	116.483805	1478794260
39.907752	116.483802	1478794266
39.908785	116.483792	1478794272
39.909615	116.483785	1478794278
39.910577	116.483807	1478794284
39.911385	116.483813	1478794290
39.912340	116.483827	1478794296
39.913268	116.483845	1478794302
39.914215	116.483843	1478794308
39.915212	116.483832	1478794314
39.916267	116.483815	1478794320
39.917197	116.483847	1478794326
39.918115	116.483827	1478794332
39.919403	116.483815	1478794338
39.920547	116.483793	1478794344
39.921457	116.483807	1478794350
39.922705	116.483822	1478794356
39.923673	116.483835	1478794362
39.925725	116.483828	1478794368
39.926870	116.483817	1478794374
39.926870	116.483817	1478794380
39.928135	116.483818	1478794386
39.930070	116.483838	1478794392
39.931182	116.483817	1478794398
39.932148	116.48381	1478794404
39.933502	116.483813	1478794410
39.934432	116.483815	1478794416
39.935387	116.483807	1478794422
39.936513	116.48379	1478794428
39.937668	116.483802	1478794434
39.938607	116.483807	1478794440
39.939718	116.483827	1478794446
39.940663	116.483827	1478794452
39.941772	116.48382	1478794458
39.943017	116.483815	1478794464
39.944130	116.483813	1478794470
39.945430	116.48381	1478794476
39.946395	116.483802	1478794482
39.946395	116.483802	1478794488
39.948492	116.483782	1478794494
39.949538	116.483755	1478794500
39.950600	116.48365	1478794506
39.951498	116.48354	1478794512
39.952388	116.483375	1478794518
39.954267	116.482807	1478794524
39.955262	116.482365	1478794530
39.956282	116.481793	1478794536
39.957282	116.481085	1478794542
39.958085	116.480387	1478794548
39.958985	116.479525	1478794554
39.959767	116.47869	1478794560
39.960395	116.47799	1478794566
39.960982	116.477225	1478794572
39.961550	116.476467	1478794578
39.962128	116.475527	1478794584
39.962665	116.474652	1478794590
39.963367	116.473545	1478794596
39.964187	116.472257	1478794602
39.964785	116.47133	1478794608
39.965618	116.470047	1478794614
39.966198	116.469163	1478794620
39.966917	116.468077	1478794626
39.967813	116.466717	1478794632
39.968628	116.465518	1478794638
39.969275	116.46456	1478794644
39.969948	116.463523	1478794650
39.970610	116.462492	1478794656
39.971263	116.461463	1478794662
39.971817	116.460595	1478794668
39.972515	116.45952	1478794674
39.973238	116.45842	1478794680
39.973852	116.457487	1478794686
39.974547	116.45645	1478794692
39.975357	116.455172	1478794698
39.976050	116.454095	1478794704
39.976637	116.453182	1478794710
39.977378	116.452075	1478794716
39.978040	116.451082	1478794722
39.978847	116.449857	1478794728
39.979537	116.448843	1478794734
39.980230	116.447808	1478794740
39.980820	116.446917	1478794746
39.981635	116.445647	1478794752
39.982337	116.44458	1478794758
39.983088	116.44342	1478794764
39.983752	116.442433	1478794770
39.984388	116.441278	1478794776
39.984792	116.440257	1478794782
39.985352	116.438748	1478794788
39.985862	116.43741	1478794794
39.986367	116.436148	1478794800
39.986917	116.43471	1478794806
39.987328	116.433567	1478794812
39.987608	116.432575	1478794818
39.987787	116.431508	1478794824
39.987853	116.430045	1478794830
39.987805	116.428603	1478794836
39.987763	116.427177	1478794842
39.987680	116.425548	1478794848
39.987603	116.424162	1478794854
39.987558	116.422782	1478794860
39.987552	116.421395	1478794866
39.987515	116.419988	1478794872
39.987465	116.418822	1478794878
39.987420	116.417158	1478794884
39.987402	116.41573	1478794890
39.987378	116.414057	1478794896
39.987337	116.412847	1478794902
39.987318	116.411362	1478794908
39.987268	116.409773	1478794914
39.987267	116.40845	1478794920
39.987255	116.407347	1478794926
39.987232	116.406135	1478794932
39.987235	116.404905	1478794938
39.987212	116.40391	1478794944
39.987270	116.40314	1478794950
39.987608	116.402677	1478794956
39.987952	116.402262	1478794962
39.988378	116.401708	1478794968
39.988825	116.401415	1478794974
39.989480	116.401417	1478794980
39.990110	116.401427	1478794986
39.990950	116.401397	1478794992
39.991850	116.401383	1478794998
39.992540	116.401367	1478795004
39.992950	116.401352	1478795010
39.993212	116.401343	1478795016
39.993212	116.401343	1478795022
39.993212	116.401343	1478795028
39.993212	116.401343	1478795034
39.993230	116.401423	1478795040
39.993323	116.401417	1478795046
39.993323	116.401417	1478795052
39.993323	116.401417	1478795058
39.993323	116.401417	1478795064
39.993323	116.401417	1478795070
39.993323	116.401417	1478795076
39.993323	116.401417	1478795082
39.993323	116.401417	1478795088
39.993323	116.401417	1478795094
39.993323	116.401417	1478795100
39.993323	116.401417	1478795106
39.993323	116.401417	1478795112
39.993417	116.401368	1478795118
39.993557	116.400988	1478795124
39.993523	116.400372	1478795130
39.993482	116.399398	1478795136
39.993447	116.398485	1478795142
39.993397	116.397705	1478795148
39.993385	116.397072	1478795154
39.993428	116.396227	1478795160
39.993368	116.395358	1478795166
39.993327	116.394415	1478795172
39.993292	116.393623	1478795178
39.993292	116.393623	1478795184
39.993292	116.393623	1478795190
39.993292	116.393623	1478795196
39.993292	116.393623	1478795202
39.993292	116.393623	1478795208
39.993292	116.393623	1478795214
39.993292	116.393623	1478795220
39.993292	116.393623	1478795226
39.993292	116.393623	1478795232
39.993292	116.393623	1478795238
39.993292	116.393623	1478795244
39.993292	116.393623	1478795250
39.993292	116.393623	1478795256
39.993292	116.393623	1478795262
39.993292	116.393623	1478795268
39.993292	116.393623	1478795274
39.993292	116.393623	1478795280
39.993292	116.393623	1478795286
39.993292	116.393623	1478795292
39.993292	116.393623	1478795298
39.993292	116.393623	1478795304
39.993292	116.393623	1478795310
39.993292	116.393623	1478795316
39.993292	116.393623	1478795322
39.993292	116.393623	1478795328
39.993292	116.393623	1478795334
39.993292	116.393623	1478795340
39.993292	116.393623	1478795346
39.993292	116.393623	1478795352
39.993292	116.393623	1478795358
39.993292	116.393623	1478795364
39.993292	116.393623	1478795370
39.993292	116.393623	1478795376
39.993292	116.393623	1478795382
39.993292	116.393623	1478795388
39.993292	116.393623	1478795394];

road_path = t.toRoads(trajectory);

% road_path has 4 column: timeSlot, timeStart, roadId, travelTime.
% timeStart is the timestamp at which the car enter the road.
% roadId is only unique for one import, if you import map a second time, it may change.
% travelTime is the time the car spent to travel through the road (seconds).

matched = t.exactPoints(trajectory);
% you got a matrix with these columns: road_id, matched_lat, matched_lon, origin_lat, origin_lon, origin_time
% you got something like this, not all points are matched.
% road_id, matched_lat, matched_lon, origin_lat, origin_lon, origin_time
% 91558	39.8621799403485	116.471331237138	39.8621800000000	116.471353000000	1478793726.00000
% 91558	39.8628920583424	116.471329285239	39.8628920000000	116.471308000000	1478793744.00000
% 91558	39.8634620705058	116.471327722848	39.8634620000000	116.471302000000	1478793756.00000
% 8733	39.8642863311037	116.471309927235	39.8642850000000	116.471250000000	1478793762.00000
% 8733	39.8647302893028	116.471300066036	39.8647280000000	116.471197000000	1478793768.00000
% 8733	39.8658407862976	116.471275399677	39.8658400000000	116.471240000000	1478793780.00000
% 8733	39.8662268850230	116.471266823652	39.8662270000000	116.471272000000	1478793792.00000
% 8734	39.8666830289700	116.471226998921	39.8666830000000	116.471233000000	1478793846.00000
% 8734	39.8672570107662	116.471229769803	39.8672570000000	116.471232000000	1478793852.00000
% 8734	39.8678978268708	116.471232863325	39.8678980000000	116.471197000000	1478793858.00000
% 8734	39.8685998829218	116.471236252482	39.8686000000000	116.471212000000	1478793864.00000
% 8734	39.8691949559484	116.471239125180	39.8691950000000	116.471230000000	1478793870.00000
% 8734	39.8697481120190	116.471241795526	39.8697480000000	116.471265000000	1478793876.00000
% 57875	39.8698396722739	116.471486794774	39.8699020000000	116.471512000000	1478793882.00000