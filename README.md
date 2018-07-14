# map-matching-matlab
simple map-matching code for GraphHopper's map-matching lib

# Usage

### STEP1
first, import map file into a dir. (You have to done this by Java directly because doing this in matlab may cause some err).
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
t.load('D:\map-matching-matlab\cache', 20) # 20 is the gps accuracy, unit one meter
# example trajectory, first column is latitude, second is longitude, third is timestamp (seconds).
trajectory=[
39.8621800000000	116.471353000000	1478793726
39.8621550000000	116.471312000000	1478793732
39.8623980000000	116.471278000000	1478793738
39.8628920000000	116.471308000000	1478793744
39.8628920000000	116.471308000000	1478793750
39.8634620000000	116.471302000000	1478793756
39.8642850000000	116.471250000000	1478793762
39.8647280000000	116.471197000000	1478793768
39.8647280000000	116.471197000000	1478793774
39.8658400000000	116.471240000000	1478793780
39.8658400000000	116.471240000000	1478793786
39.8662270000000	116.471272000000	1478793792
39.8662270000000	116.471272000000	1478793798
39.8662270000000	116.471272000000	1478793804
39.8662270000000	116.471272000000	1478793810
39.8662270000000	116.471272000000	1478793816
39.8662270000000	116.471272000000	1478793822
39.8662270000000	116.471272000000	1478793828
39.8662270000000	116.471272000000	1478793834
39.8662270000000	116.471272000000	1478793840
39.8666830000000	116.471233000000	1478793846
39.8672570000000	116.471232000000	1478793852
39.8678980000000	116.471197000000	1478793858
39.8686000000000	116.471212000000	1478793864
39.8691950000000	116.471230000000	1478793870
39.8697480000000	116.471265000000	1478793876
39.8699020000000	116.471512000000	1478793882]

road_path = t.toRoads(trajectory)
```
