map-matching-matlab
-----------------------
Use GraphHopper's map-matching lib from matlab (a java wrapper)

# Usage

### STEP1
first, compile this project and import map file into a dir. (You have to done this by Java directly because doing this in matlab may cause some err).
```
java -cp "*" edu.buaa.ACTMapMatching D:\map-matching-matlab\beijing_around.pbf D:\map-matching-matlab\cache
```

### STEP2
then, modify matlab classpath file to allow seek of the map-matching lib. Type `edit classpath.txt` in matlab console,
paste following to the file:
```
# D:/map-matching-matlab/out/artifacts/map_matching4matlab_jar/commons-compress-1.10.jar
# D:/map-matching-matlab/out/artifacts/map_matching4matlab_jar/commons-io-1.3.1.jar
# D:/map-matching-matlab/out/artifacts/map_matching4matlab_jar/commons-logging-1.0.4.jar
D:/map-matching-matlab/out/artifacts/map_matching4matlab_jar/graphhopper-core-0.8.2.jar
D:/map-matching-matlab/out/artifacts/map_matching4matlab_jar/graphhopper-map-matching-core-0.8.2.jar
D:/map-matching-matlab/out/artifacts/map_matching4matlab_jar/graphhopper-reader-osm-0.8.2.jar
D:/map-matching-matlab/out/artifacts/map_matching4matlab_jar/graphhopper-tools-lgpl-0.8.2.jar
D:/map-matching-matlab/out/artifacts/map_matching4matlab_jar/hmm-lib-1.0.0.jar
D:/map-matching-matlab/out/artifacts/map_matching4matlab_jar/map-matching4matlab.jar
D:/map-matching-matlab/out/artifacts/map_matching4matlab_jar/osmosis-osm-binary-0.44.1.jar
D:/map-matching-matlab/out/artifacts/map_matching4matlab_jar/protobuf-java-2.6.1.jar
# D:/map-matching-matlab/out/artifacts/map_matching4matlab_jar/slf4j-api-1.7.7.jar
D:/map-matching-matlab/out/artifacts/map_matching4matlab_jar/trove4j-3.0.3.jar
# D:/map-matching-matlab/out/artifacts/map_matching4matlab_jar/xmlgraphics-commons-2.1.jar#
# D:\matlabcalljava\out\artifacts\matlabcalljava_jar\matlabcalljava.jar
```
some are commented because matlab already have include such jar (although it could be a different version).
restart matlab to make above config work.

### STEP3
finally, call java methods from matlab. (check `example.m` for example)
```
t = edu.buaa.ACTMapMatching
t.load('D:\map-matching-matlab\cache', 20); % 20 is the gps accuracy, unit one meter
% example trajectory, first column is latitude, second is longitude, third is timestamp (seconds).
trajectory=[
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
39.869902	116.471512	1478793882];

road_path = t.toRoads(trajectory);
% you got something like this,
%
% 0	1478793829.489	8733	73
% 0	1478793878.839	8734	49
%
% the second column is the timestamp (second) at which the car enter the road, the first column is the timeslot the timestamp falls (48 timeslot a day)
% the third column is roadId which is only unique for one import, if you import map a second time, it may change.
% the last column is travelTime which is the time the car spent to travel through the road (seconds).


matched = t.exactPoints(trajectory);
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
```
