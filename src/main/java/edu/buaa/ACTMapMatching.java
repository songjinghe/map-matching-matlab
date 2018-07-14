package edu.buaa;

import com.graphhopper.matching.EdgeMatch;
import com.graphhopper.matching.GPXExtension;
import com.graphhopper.matching.MapMatching;
import com.graphhopper.matching.MatchResult;
import com.graphhopper.reader.osm.GraphHopperOSM;
import com.graphhopper.routing.AlgorithmOptions;
import com.graphhopper.routing.util.CarFlagEncoder;
import com.graphhopper.routing.util.EncodingManager;
import com.graphhopper.routing.weighting.FastestWeighting;
import com.graphhopper.routing.weighting.Weighting;
import com.graphhopper.util.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;
import java.util.TimeZone;

public class ACTMapMatching {

    public static DistanceCalc distanceCalc = new DistancePlaneProjection();
    private GraphHopperOSM hopper = new GraphHopperOSM();
    public MapMatching mapMatching;
    private File cacheDir;

    public static void main(String[] args)
    {
        if(args.length>1){
            System.out.println("import map data from "+args[0]+" to "+args[1]);
            importMap(args[0], args[1]);
        }else{
            System.err.println("args too less: need map_file_path and cache_dir_path");
        }
    }

    private static void importMap(String mapFileStr, String cacheDirStr){
        File mapFile = new File(mapFileStr);
        File cacheDir = new File(cacheDirStr);
        if(!mapFile.exists()){
            System.err.println("map file not exist!");
            return;
        }
        if(!cacheDir.exists()){
            try {
                Files.createDirectories(cacheDir.toPath());
            } catch (IOException e) {
                System.err.println(e.getMessage());
                return;
            }
        }else if(cacheDir.isFile()){
            System.err.println("invalid cache dir! need folder but find file.");
            return;
        }
        GraphHopperOSM hopper = new GraphHopperOSM();
        CarFlagEncoder encoder = new CarFlagEncoder();
        // import OpenStreetMap data
        hopper.setDataReaderFile(mapFile.getAbsolutePath());
        hopper.setGraphHopperLocation(cacheDir.getAbsolutePath());
        hopper.setEncodingManager(new EncodingManager(encoder));
        hopper.getCHFactoryDecorator().setEnabled(false);
        hopper.importOrLoad();
    }

    public void load(String tmpDir, double gpsAccuracy)
    {
        File cacheDir = new File(tmpDir);
        if(!cacheDir.exists()){
            System.err.println("invalid cache dir! folder not exist.");
        }else if(cacheDir.isFile()){
            System.err.println("invalid cache dir! need folder but find file.");
            return;
        }else if(cacheDir.list()==null){

        }
        this.cacheDir = cacheDir;
        CarFlagEncoder encoder = new CarFlagEncoder();
        // import OpenStreetMap data
//        hopper.setDataReaderFile(mapFile.getAbsolutePath());
//        hopper.setGraphHopperLocation(cacheDir.getAbsolutePath());
        hopper.setEncodingManager(new EncodingManager(encoder));
        hopper.getCHFactoryDecorator().setEnabled(false);
        hopper.load(cacheDir.getAbsolutePath());
        // create MapMatching object, can and should be shared across threads
        String algorithm = Parameters.Algorithms.DIJKSTRA_BI;
        Weighting weighting = new FastestWeighting(encoder);
        AlgorithmOptions algoOptions = AlgorithmOptions.start()
                .maxVisitedNodes(2100000000)
                .algorithm(algorithm)
                .weighting(weighting)
                .build();

        mapMatching = new MapMatching(hopper, algoOptions);
        mapMatching.setMeasurementErrorSigma(gpsAccuracy);
    }

    public void clearMap(){
        try {
            delete(this.cacheDir);
            System.out.println("map cache cleaned.");
        }catch (IOException e){
            System.err.println("clean failed, try delete '"+this.cacheDir.getAbsolutePath()+"' manually.");
        }
    }

    private void delete(File f) throws IOException {
        if (f.isDirectory()) {
            for (File c : f.listFiles()) delete(c);
        }
        if (!f.delete()) throw new FileNotFoundException("Failed to delete file: " + f);
    }

    public double[][] toRoads(double[][] traj) {
        List<GPXEntry> data = input(traj);
        try {
            return matching(data);
        } catch (IllegalStateException e) {
            System.err.println("MM failed: Unknown error in map-matching lib.");
        } catch (NoEdgeMatchedException e) {
            System.err.println("MM failed: No road matched.");
        } catch (CannotCalcTravelTimeException e) {
            System.err.println("MM failed: Unable to calculate time.");
        } catch (IllegalArgumentException e) {
            System.err.println("MM failed: Seems get lost.");
        } catch (RuntimeException e) {
            if (e.getMessage() != null && e.getMessage().startsWith("Sequence is broken for submitted track at time step")) {
                System.err.println("MM failed: Too long to match.");
            }else{
                System.err.println("MM failed: Runtime error in map-matching lib.");
            }
        }
        return null;
    }

    public double[][] exactPoints(double[][] traj) {
        List<GPXEntry> data = input(traj);
        try {
            MatchResult mr = mapMatching.doWork(data);
            List<EdgeMatch> matches = mr.getEdgeMatches();
            if(matches.isEmpty()){
                throw new NoEdgeMatchedException();
            }else {
                List<double[]> result = new ArrayList<>();
                for (EdgeMatch edge : matches) {
                    int edgeID = edge.getEdgeState().getEdge();
                    List<GPXExtension> gpsPointOnTheWay = edge.getGpxExtensions();
                    int k = 0;
                    PointList nodes = edge.getEdgeState().fetchWayGeometry(3);
                    for (int j = 0; j < nodes.size() - 1; j++) { // loop through edges. edgeNumber = nodeNumber - 1
                        double startLat = nodes.getLatitude(j);
                        double startLon = nodes.getLongitude(j);
                        double endLat = nodes.getLatitude(j + 1);
                        double endLon = nodes.getLongitude(j + 1);
                        ProjectionResult r;
                        while (k < gpsPointOnTheWay.size()) {
                            GPXEntry gpxPoint = gpsPointOnTheWay.get(k).getEntry();
                            r = projectionPoint(
                                    startLat, startLon,
                                    endLat, endLon,
                                    gpxPoint.getLat(), gpxPoint.getLon());
                            result.add( new double[]{edgeID, r.y, r.x, gpxPoint.getLat(), gpxPoint.getLon(), gpxPoint.getTime()/1000} );
                        }
                    }
                }
                double[][] toReturn = new double[result.size()][];
                int i=0;
                for(double[] line : result){
                    toReturn[i] = line;
                    i++;
                }
                return toReturn;
            }
        } catch (IllegalStateException e) {
            System.err.println("MM failed: Unknown error in map-matching lib.");
        } catch (NoEdgeMatchedException e) {
            System.err.println("MM failed: No road matched.");
        } catch (CannotCalcTravelTimeException e) {
            System.err.println("MM failed: Unable to calculate time.");
        } catch (IllegalArgumentException e) {
            System.err.println("MM failed: Seems get lost.");
        } catch (RuntimeException e) {
            if (e.getMessage() != null && e.getMessage().startsWith("Sequence is broken for submitted track at time step")) {
                System.err.println("MM failed: Too long to match.");
            }else{
                System.err.println("MM failed: Runtime error in map-matching lib.");
            }
        }
        return null;
    }

    private List<GPXEntry> input(double[][] traj){
        List<GPXEntry> result = new ArrayList<>();
        for(int i=0; i<traj.length; i++){
            double[] row = traj[i];
            result.add(new GPXEntry(row[0], row[1], ((long) row[2])*1000));
        }
        return result;
    }

    private double[][] matching(List<GPXEntry> traj) throws NoEdgeMatchedException {
        MatchResult mr = mapMatching.doWork(traj);

        List<EdgeMatch> matches = mr.getEdgeMatches();
        if(matches.isEmpty()){
            throw new NoEdgeMatchedException();
        }else {
            List<TrajectoryRoadEntry> edges = calcTravelTime(matches);
            return output(edges);
        }
    }

    private double[][] output(List<TrajectoryRoadEntry> edges) {
        double[][] result = new double[edges.size()][];
        int i=0;
        for(TrajectoryRoadEntry e : edges){
            result[i] = new double[]{e.timeSlot, e.timeStart/1000, e.edgeID, e.travelTime};
            i++;
        }
        return result;
    }

    /**
     * calculate travel time of each road along the trajectory.
     * we need this function because the HMM map-matching algorithm (used by graphHopper)
     * do not generate the travel time of each road.
     *
     * First, we transform each known gps entry point in the trajectory from [latitude, longitude, timestamp] to
     * [distance, timestamp], where the 'distance' is the distance (not Displacement) from current gps point
     * to the trajectory start gps point. We also transform query points (joint point of road) to distance.
     * Then we use linear interpolation to calculate the timestamp of query point.
     * Finally, the difference between the two timestamp of the joint point of road is the travel time.
     *
     *
     * @param matches map-matching result
     * @return [roadID, timeSlot, travelTime] list
     */
    private List<TrajectoryRoadEntry> calcTravelTime(List<EdgeMatch> matches)
    {
        List<DistanceTimeEntry> knownPoints = new ArrayList<DistanceTimeEntry>();
        List<DistanceTimeEntry> queryPoints = new ArrayList<DistanceTimeEntry>();

        // collection all available [distance v.s. timestamp] data.
        double trajLen = calcDistanceTimeList(matches, knownPoints, queryPoints);

//        traj.setLength(trajLen);

        if(knownPoints.size()<2) throw new CannotCalcTravelTimeException();

        // calc joint point timestamp from its distance(to trajectory start point) using linear interpolation.
        try {
            return calcTravelTime(knownPoints, queryPoints);
        }catch (CannotCalcTravelTimeException e){
//            long id = System.currentTimeMillis();
//            traj.toFile(new File(id+".gpx"));
            throw e;
        }
    }

    // @return traj length (in meter)
    private double calcDistanceTimeList(List<EdgeMatch> matches, List<DistanceTimeEntry> knownPoints, List<DistanceTimeEntry> queryPoints)
    {
        double disToTrajectoryStart = 0;
        for (EdgeMatch edge : matches) {
            List<GPXExtension> gpsPointOnTheWay = edge.getGpxExtensions();
            int k = 0;
            double disToWayStart = 0;
            PointList nodes = edge.getEdgeState().fetchWayGeometry(3);
            for (int j = 0; j < nodes.size() - 1; j++) { // loop through edges. edgeNumber = nodeNumber - 1
                double startLat = nodes.getLatitude(j);
                double startLon = nodes.getLongitude(j);
                double endLat = nodes.getLatitude(j + 1);
                double endLon = nodes.getLongitude(j + 1);
                ProjectionResult r;
                while (k < gpsPointOnTheWay.size()) {
                    GPXEntry gpxPoint = gpsPointOnTheWay.get(k).getEntry();
                    r = projectionPoint(
                            startLat, startLon,
                            endLat, endLon,
                            gpxPoint.getLat(), gpxPoint.getLon());
                    if (r.inside) {
                        double deltaDis = distanceCalc.calcDist(r.y, r.x, startLat, startLon);
                        double totalDis = deltaDis + disToWayStart + disToTrajectoryStart;
                        knownPoints.add(new DistanceTimeEntry(totalDis, gpxPoint.getTime()));
                        k++;
//                        System.out.printf("%f,%f\n", r.y, r.x);
//                        System.out.printf("%f,%d,know\n", totalDis, gpxPoint.getTime());
                    } else {
                        break;
                    }
                }
                double edgeLength = distanceCalc.calcDist(startLat, startLon, endLat, endLon);
                disToWayStart += edgeLength;
//                // query all nodes (include pillar node)
//                queryPoints.add(disToTrajectoryStart+disToWayStart);
//                queryIndex.add(knownPoints.size());
            }
            disToTrajectoryStart += disToWayStart;
            // only query tower nodes.
            queryPoints.add(new DistanceTimeEntry(disToTrajectoryStart, edge.getEdgeState().getEdge()));
//            queryIndex.add(knownPoints.size());
        }
        return disToTrajectoryStart;
    }

    // calc joint point timestamp from its distance(to trajectory start point) using linear interpolation.
    private List<TrajectoryRoadEntry> calcTravelTime(List<DistanceTimeEntry> known, List<DistanceTimeEntry> query)
    {
        List<TrajectoryRoadEntry> results = new ArrayList<TrajectoryRoadEntry>();
        int k=0;
        for(int i=0; i<query.size(); i++)
        {
            DistanceTimeEntry cur = query.get(i);
            while(k < known.size() && known.get(k).distance < cur.distance){ k++; }

            if(0<k && k<known.size()) // normal
            {
                DistanceTimeEntry pre = known.get(k-1);
                DistanceTimeEntry post = known.get(k);
                cur.time = (long) ( pre.time +
                        (post.time - pre.time) *
                                (cur.distance - pre.distance)/
                                (post.distance - pre.distance));
                if(i>0)
                {
                    DistanceTimeEntry last = query.get(i-1);
                    if(last.time!=0) {
                        long travelTime = (cur.time - last.time) / 1000;
                        if (travelTime < 0 ) {
                            System.err.printf("invalid travel time(%d)! cur: %d pre: %d\n", travelTime, cur.time, last.time);
                        }else if( travelTime > 1800 ){
                            throw new CannotCalcTravelTimeException();
                        }else {
                            results.add(new TrajectoryRoadEntry(
                                    cur.edgeID,
                                    cur.time,
                                    timestamp2TimeSlot(cur.time),
                                    (int) travelTime));
                        }
                    }
                }
            }else if(k==known.size()){
                // query a point after(k==know.size()) all know points. exit loop.
                break;
            }else{
                // query a point before(k==0) all know points. do nothing
            }
        }
        return results;
    }

    private int timestamp2TimeSlot(long timestamp)
    {
        Calendar c = Calendar.getInstance(TimeZone.getTimeZone("GMT+8:00"));
        c.setTimeInMillis(timestamp);
        int timeSlot = c.get(Calendar.HOUR_OF_DAY)*2;
        if(c.get(Calendar.MINUTE)>=30){
            timeSlot++;
        }
        return timeSlot;
    }

    static ProjectionResult projectionPoint(double startY, double startX, double endY, double endX, double y, double x)
    {
//            System.out.printf("(%f %f) (%f %f) (%f %f) \n", startY, startX, endY, endX, y, x);
        double roadX = endX - startX;
        double roadY = endY - startY;
        double len = Math.hypot(roadX, roadY);
        double iX = roadX / len;
        double iY = roadY / len;
        double vX = x-startX;
        double vY = y-startY;
        double innerProduct = vX * roadX + vY * roadY;
        double projectionLen = innerProduct / len;
        double tX = startX + projectionLen * iX;
        double tY = startY + projectionLen * iY;
        return new ProjectionResult(tX, tY, projectionLen<=len);
    }

    private static class DistanceTimeEntry
    {
        public double distance;
        public long time;
        public int edgeID=-1;
        public DistanceTimeEntry(double dis, long t)
        {
            this.distance = dis;
            this.time = t;
        }
        public DistanceTimeEntry(double dis, int edgeID)
        {
            this.distance = dis;
            this.edgeID = edgeID;
        }
    }

    static class ProjectionResult
    {
        public double x;
        public double y;
        public boolean inside;
        public ProjectionResult(double x, double y, boolean b)
        {
            this.x = x;
            this.y = y;
            this.inside = b;
        }
    }

    public static class TrajectoryRoadEntry
    {
        public final int edgeID;
        public final long timeStart;
        public final int timeSlot;
        public final int travelTime; // in seconds;

        public TrajectoryRoadEntry(int edgeID, long timeStart, int timeSlot, int travelTime)
        {
            this.edgeID = edgeID;
            this.timeStart = timeStart;
            this.timeSlot = timeSlot;
            this.travelTime = travelTime;
        }

        public String toString(){
            return String.format("%d %d %d", edgeID, timeStart, travelTime);
        }
    }


    public void help()
    {
        System.out.println("");
    }

    private class CannotCalcTravelTimeException extends RuntimeException
    {
    }

    private class NoEdgeMatchedException extends Throwable {
    }
}


