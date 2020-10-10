package edu.buaa;

import com.graphhopper.util.GPXEntry;

import java.util.*;
import java.util.concurrent.Semaphore;

public class KdAnonymous
{
    private static void printf(String str) {
        if(str.length()>0 && str.charAt(str.length()-1)=='\n') str = str.substring(0, str.length()-1);
        System.err.println(str);
    }
    private static double SQUARE(double x){
        return x*x;
    }
    private static void exit(int code) {
        throw new RuntimeException();
    }


    // nD-Objects.h
    private static final int MAX_POINTS_PER_OBJECT = 50000000;
    private static class database_nD {
        int n;
        int dim;
        double[][] x;
        int[] name;
    }
    private static class txy_type{
        int t;
        double x,y;
    }
    private static class trajectory{
        int name;
        int from, to; /* start and end time */
        double[] x; /* nD-coordinates (n=2*(to-from+1)) */
        public String toString(){
            return name+" "+from+" "+to+" "+x.length;
        }
    }
    private static class database_trajectory{
        int n; // number of trajectories.
        int allocated_space;
        trajectory[] db;
        double pseudo_diameter; /* Estimate of dataset diameter (Actual = diagonal of the MBB)*/
        int max_eq_class_size;
        int max_traj_dim;
        int traj_deleted;
        int n_points;
        int points_deleted; /* Including whole deleted trajectories */
    }

    //database_trajectory load_objects_compact(char * nomefile, int break_every_n);
    private static database_trajectory load_objects_compact(txy_type[] txy_buffer, List<List<GPSPoint>> data, int break_every_n){
        boolean go_on, process_object;
        int id, prev_id;
        int buffer_size;
        int t;
        double x,y;
        int i, new_i, j, k;
        //txy_type txy_buffer[MAX_POINTS_PER_OBJECT];
        txy_type txy_tmp;
        double[] new_entry;
        double xmin=Double.MAX_VALUE, xmax=Double.MIN_VALUE, ymin=Double.MAX_VALUE, ymax=Double.MIN_VALUE; /* MBB of dataset*/
        database_trajectory output = new database_trajectory();// = { 0, 0, null, 0, 0, 0, 0, 0, 0 };


        printf("  Loading objects...");
//        ff=fopen(nomefile, "r");
//        if(ff == NULL) { printf("Error: cannot open input file!\n"); exit(-1);};

        /* Init output database */
        output.n=0;
        output.allocated_space=data.size();
        output.db = new trajectory[output.allocated_space];
        /* flags */
        go_on=true;
        process_object=false;
        prev_id=-999999;
        buffer_size=0;


        for(int index=0;index<data.size();index++){
            List<GPSPoint> list = data.get(index);
            trajectory traj = new trajectory();
            traj.name = list.get(0).getCarID();
            traj.from = list.get(0).getTimeSlot();
            traj.to = list.get(list.size()-1).getTimeSlot();
            traj.x = new double[list.size()*2];
            for(int iii=0; iii<list.size(); iii++){
                GPSPoint p = list.get(iii);
                traj.x[2*iii] = p.getLat();
                traj.x[2*iii+1] = p.getLon();
                xmin=Math.min(xmin, p.getLat());
                xmax=Math.max(xmax, p.getLat());
                ymin=Math.min(ymin, p.getLon());
                ymax=Math.max(ymax, p.getLon());
            }
            output.db[index] = traj;
        }
        output.n = data.size();
        output.traj_deleted = 0;
        output.points_deleted = 0;



//        while (go_on)
//        {
//            /*  Reads a line of original Brinhoff's format: */
//            if(fscanf(ff, "%d\t%d\t%lf\t%lf",
//                    &id,&t,&x,&y)<4)
//            { /* EOF or truncated line  =>  stop */
//                go_on=false; /* stop reading file */
//                process_object=true; /* process the data collected so far */
//            }
//            else
//            {
//                /* Line was read correctly */
//                if(output.n_points==0)
//                    prev_id=id;
//                if(id!=prev_id) /* end-of-previous object */
//                    process_object=true;
//                output.n_points++;
//
//                /* Stats for computing pseudo_diameter */
//                xmin=Math.min(xmin,x);
//                xmax=Math.max(xmax,x);
//                ymin=Math.min(ymin,y);
//                ymax=Math.max(ymax,y);
//            }
//
//            if(process_object)
//            { /* The buffer content has to be converted in a trajectory */
//
//                /* Step 1: sort the data w.r.t. time */
//                /* TODO: replace bubble sort with qsort or something else */
//                for(i=buffer_size-1; i>0;)
//                {
//                    new_i=0;
//                    for(j=0; j<i; j++)
//                        if(txy_buffer[j].t>txy_buffer[j+1].t)
//                        {
//                            txy_tmp=txy_buffer[j];
//                            txy_buffer[j]=txy_buffer[j+1];
//                            txy_buffer[j+1]=txy_tmp;
//                            new_i=j;
//                        }
//                    i=new_i;
//                }
//                /* Step 2: cut head and tail, if needed */
//                for(i=0; (i<buffer_size) && ((txy_buffer[i].t % break_every_n) != 0); i++);
//                for(j=buffer_size-1; (j>0) && ((txy_buffer[j].t % break_every_n) != 0); j--);
//                if(i<=j) /* at least one point was not cut => process points in [i,j] */
//                {
//                    /* Step 3: make it a 2N-dim vector for a nD_database */
//                    new_entry= new double[2*(j-i+1)];
////                    if(new_entry==NULL)
////                    {
////                        printf("call `malloc(%d)` returns NULL! @_@!!", 2*(j-i+1)*sizeof(double));
////                        exit(-250);
////                    }
//                    for(k=i; k<=j; k++)
//                    {
//                        new_entry[2*(k-i)]=txy_buffer[k].x;
//                        new_entry[2*(k-i)+1]=txy_buffer[k].y;
//                    }
//                    /* Step 4: add the new entry to trajectory repository */
//                    if(output.n>=output.allocated_space) /* realloc space for "output" if needed */
//                    {
//                        if(output.n>output.allocated_space) { printf("ERROR: n>allocated_space!\n"); exit(-1);}
//                        k=(int)(output.allocated_space * 1.2);
//                        output.db= new trajectory[k];
//                        output.allocated_space=k;
//                    }
//                    output.db[output.n].x=new_entry;
//                    output.db[output.n].from=txy_buffer[i].t;
//                    output.db[output.n].to=txy_buffer[j].t;
//                    output.db[output.n].name=prev_id;
//                    output.n++;
//                    output.points_deleted+=buffer_size-(j-i+1);
//                } else {
//                    output.traj_deleted++;
//                    output.points_deleted+=buffer_size; /* NOTE: points_deleted includes whole trajectories */
//                }
//                buffer_size=0; /* Empty buffer */
//                process_object=false;
//                prev_id=id;
//            }
//            if(go_on)
//            {
//                /* Put the data read (which had no errors) in the buffer */
//                if(buffer_size>=MAX_POINTS_PER_OBJECT)
//                { printf(String.format("ERROR: MAX_POINTS_PER_OBJECT (%d) exceeded !\n", MAX_POINTS_PER_OBJECT)); exit(-1); }
//                txy_buffer[buffer_size].t=t;
//                txy_buffer[buffer_size].x=x;
//                txy_buffer[buffer_size].y=y;
//                buffer_size++;
//            }
//        };

        output.pseudo_diameter = Math.sqrt(SQUARE(xmax-xmin)+SQUARE(ymax-ymin));
        /*
          printf("###############\nxmin %lf  xmax %lf  ymin %lf  ymax %lf  pseudo-diam %lf\n###############\n",xmin,xmax,ymin,ymax,output.pseudo_diameter);
        */

        printf(" Done.\n");
        if(output.n>0) {
            printf(String.format("  -> Trajectories: %d, Points: %d, Diameter: %.3f", output.n, output.n_points, output.pseudo_diameter));
            printf(String.format("  -> Removed Trajectories: %d, Removed Points: %d", output.traj_deleted, output.points_deleted));
        }
        return output;
    }

    private static FasterEuclidean_nD disFaster = null;
    private static double euclidean_nD(database_nD DB, int a, int b){
        return original_euclidean_nD(DB, a, b);
    }

    private static double original_euclidean_nD(database_nD DB, int a, int b){
        double x1,x2,sum;
        int i;

        for(i=0, sum=0; i<DB.dim; i++)
        {
            x1=DB.x[a][i]; /* No range check -> faster and allows use of dummy objects */
            x2=DB.x[b][i];
            sum+=(x1-x2)*(x1-x2);
        }
        //        System.out.println(a+" "+b+" "+dis);
        return Math.sqrt(sum);
    }

    /**
     *  nwa.h
     */

    private static class heap_obj {
        int i;
        double d;

        public heap_obj(int i, double d) {
            this.i = i;
            this.d = d;
        }
    }

    private static class output_clusters_struct {
        int n;
        int[]  c_size;
        int[][] c_objects;
    }

    // union
    private static class foint implements Comparable<foint>{
        heap_obj v; // void *

        public foint(int j, double v) {
            this.v = new heap_obj(j,v);
        }

        @Override
        public int compareTo(foint o) {
            return Double.valueOf(o.v.d).compareTo(this.v.d);
        }
    } ;


//    /**
//     *  NEVER_WALK_ALONE clustering algorithm ORIGINAL
//     */
//    private static output_clusters_struct never_walk_alone(database_nD DB, int K, double first_max_radius, int first_max_trash_size)
//    {
//        /*  NOTICE: DB->x has to contain at least DB.n+1 vectors -- the extra one is a "dummy trajectory"
//            The (DB.n+1)-th should be allocated as the others, though with any values */
//
//        PriorityQueue<foint> myHeap; //smallest at top;
//        foint myfoint;
//        int i, S_size, pivot, best_pivot, j, last_cluster = 0, next_i;
//        int[] assignments;
//        int[] cluster_size;
//        int[] alive;
//        int[] pivots;
//        byte[] ex_pivots;
//        int n_ex_pivots;
//        int[] k_buffer;
//        double best_d, next_d, tmp_d, max_radius, farthest_d;
//        int max_trash_size;
//        double[] center;
//        output_clusters_struct output_clusters = new output_clusters_struct();
//        boolean solution_found;
//
//        if(DB.n<K)
//        {
//            printf(String.format("ERROR: Never_walk_alone called with K (%d) greater than database_size (%d)!\n", K, DB.n));
//            exit(-1);
//        }
//
//        /* Global Memory allocations and Initializations */
//
//        /* Objects can be: (1) alive; (2) alive but already used as pivots - and failed; (3) dead */
//        alive = new int[DB.n];       /* alive[0..S_size-1] contains IDs of alive objects */
//        ex_pivots = new byte[DB.n]; /* ex_pivots[i]==0 <=> object "i" was never used as pivot */
//        pivots = new int[DB.n];      /* pivots[i] = pivot/center of cluster "i" */
//        k_buffer = new int[K];
//        assignments = new int[DB.n]; /* assignments[i] = cluster of object "i" */
//        cluster_size = new int[DB.n]; /* (oversized, to avoid problems...) */
//
//        myHeap = new PriorityQueue<>();//HeapAlloc(K,heap_cmp_reverse); /* Heap used for efficiently collecting (K-1)-closest objects */
//
//        /* First pivot = center of the dataset --> trick: uses a dummy object of the dataset */
//        center=DB.x[DB.n];
//        for(j=0; j<DB.dim; j++)
//            center[j]=0.0;
//        for(i=0; i<DB.n; i++)
//            for(j=0; j<DB.dim; j++)
//                center[j]+=DB.x[i][j];
//        for(j=0; j<DB.dim; j++)
//            center[j] = center[j] / DB.n;
//
//        /* Cluster constraints - starting values*/
//        max_radius=first_max_radius*Math.sqrt((double)DB.dim/2.0); /* Takes into account dimensionality */
//        max_trash_size=first_max_trash_size;
//
//        /* Here the external cycle starts */
//        solution_found=false;
//        while (!solution_found)
//        {
//            /* printf(" [max_radius=%f] ",max_radius); fflush(NULL); */
//            /* Init variables */
//            for(i=0; i<DB.n;  i++)
//            {
//                alive[i]=i;
//                ex_pivots[i]=0;
//                cluster_size[i]=0; /* CHECK: not needed */
//            }
//            S_size=DB.n;   /* Number of alive objects */
//            n_ex_pivots=0;  /* Number of alive yet ex-pivot objects */
//            last_cluster=-1;
//            next_i=-1;	/* Farthest pivot from actual one (-1 = unknown) */
//            next_d=-1;
//            pivot=DB.n;	/* First pivot = center of the dataset */
//
//            /* Step 1: find clusters of size K and radius <= max_radius */
//            while ((S_size >= K) && (S_size > n_ex_pivots))
//            { /* Stop if not enough objects alive or all are ex-pivots */
//
//                last_cluster++;
//                cluster_size[last_cluster]=0;
//
//                /* Choose next pivot -- excluding ex-pivots, that cannot be pivot again */
//                if(next_i != -1)
//                {
//                    i=next_i;
//                    best_d=next_d;
//                } else
//                {
//                    best_d=-1;
//                    i=-1;
//                    for(j=0; j<S_size; j++) /* scan all alive objects */
//                        if( ex_pivots[alive[j]]==0 && ((tmp_d=euclidean_nD(DB, pivot, alive[j])) > best_d) )
//                        {
//                            i=j;
//                            best_d = tmp_d;
//                        }
//                }
//                if(i==-1) { printf("ERROR: inconsistency looking for next pivot!\n"); exit(-1); }
//                pivot=alive[i]; /* Careful: "i" has to be keeped till the end of the cycle! */
//
//                /* Select closest K-1 objects and find next possible pivot */
//                next_d=-1;
//                next_i=-1;
//                for(j=0; j<S_size; j++)
//                {
//                    if(alive[j]!=pivot)
//                    {
//                        myHeap.add(new foint(j,euclidean_nD(DB, pivot, alive[j])));
//                        if(myHeap.size()>=K) /* keep only top K-1 in the heap */
//                        {
//                            myfoint=myHeap.poll();
//                            /* Check whether the object can be next pivot */
//                            if( ( myfoint.v.d > next_d ) && ( ex_pivots[alive[myfoint.v.i]]==0 ) )
//                            {
//                                next_i=myfoint.v.i;
//                                next_d=myfoint.v.d;
//                            }
////                            free(myfoint.v);    /* Remove farthest object */
//                        }
//                    }
//                }
//
//                /* Get the distance of farthest object -- it's on the top of the heap */
//                if(myHeap.size()>0)
//                {
//                    myfoint=myHeap.peek();
//                    farthest_d=myfoint.v.d;
//                }
//                else
//                    farthest_d=0; /* Very strange case: occurs only when S_size = K = 1 */
//
//                /* Pop all the content of the heap, and copy the "i"s (indexes of alive[]) to k_buffer */
//                for(j=0; j<K-1; j++)
//                {
//                    myfoint=myHeap.poll();
//                    k_buffer[j]=myfoint.v.i;
//                }
//
//                /* Check cluster diameter (more exactly: max distance from cluster pivot) */
//                if(farthest_d <= max_radius) /* Check cluster radius (="pseudo-radius") */
//                {
//                    /* 1. "pivot" becomes a real pivot */
//                    pivots[last_cluster]=pivot;
//
//                    /* 2. From now on, check out all the members of the cluster */
//                    /* Add pivot to k_buffer, whose elements are going to be "checked out" */
//                    k_buffer[K-1]=i; /* "i" is the position of the pivot we are deleting */
//
//                    /* To avoid problems with alive[], remove objects from right to left */
//                    /* OPTIMIZE: can we avoid sorting? */
//                    Arrays.sort(k_buffer);
//                    for(j=K-1; j>=0; j--)
//                    {
//                        assignments[alive[k_buffer[j]]] = last_cluster;
//                        cluster_size[last_cluster]++;
//                        if(ex_pivots[alive[k_buffer[j]]]==1) /* remove dead ex-pivots */
//                            n_ex_pivots--;
//                        /* Virtually pop out the object from alive[] */
//                        alive[k_buffer[j]] = alive[S_size-1];
//                        S_size--;
//                        if(next_i == S_size) /* I moved my next pivot... */
//                            next_i=k_buffer[j];
//                    }
//                }
//                else
//                {
//                    /* do nothing, just "deactivate" the pivot */
//                    ex_pivots[pivot]=1;
//                    n_ex_pivots++;
//                    last_cluster--; /* undo last step */
//                }
//            }
//
//            /* Step 2: assign unassigned objects to existing clusters (where possible) */
//            /*         alive[] contains the remaining objects: try to assign them  */
//            for(i=0; i<S_size; i++)
//            {
//                if(last_cluster>=0) /* If there is at least one cluster, find closest one */
//                {
//                    best_pivot = 0;
//                    best_d=euclidean_nD(DB, pivots[0], alive[i]);
//                    for(j=1; j<=last_cluster; j++)
//                    {
//                        if((tmp_d=euclidean_nD(DB, pivots[j], alive[i])) < best_d)
//                        { best_pivot=j; best_d=tmp_d; }
//                    }
//                    if(best_d > max_radius)
//                        best_pivot = last_cluster+1; /* This is the trash */
//                }
//                else
//                    best_pivot = last_cluster+1; /* This is the trash */
//
//                assignments[alive[i]] = best_pivot;
//                cluster_size[best_pivot]++;
//            }
//
//            /* Check trash size */
//            solution_found = (cluster_size[last_cluster+1] <= max_trash_size);
//            if(!solution_found)
//            {
//                max_radius = max_radius * 1.5;		/* Exponential-growth radius threshold */
//                max_trash_size = max_trash_size;	/* Constant-valued trash threshold */
//                printf(String.format("solution not found, increase max_radius to %f", max_radius));
//            }
//        }
//
//        /* Package results for output */
//        output_clusters.n = last_cluster+2;
//        output_clusters.c_size = cluster_size;
//        output_clusters.c_objects = new int[output_clusters.n][];
//        for(i=0; i<output_clusters.n; i++)
//            output_clusters.c_objects[i] = new int[cluster_size[i]];
//        for(i=0; i<output_clusters.n; i++)
//            cluster_size[i]=0;
//        for(i=0; i<DB.n; i++)
//        {
//            output_clusters.c_objects[assignments[i]][cluster_size[assignments[i]]] = i;
//            cluster_size[assignments[i]]++;
//        }
//        return output_clusters;
//    }

    public static List<List<GPSPoint>> run(List<List<GPSPoint>> data, int anonymousK) {
        List<List<GPSPoint>> result = new ArrayList<>();

        txy_type[] txy_buffer=null;// = new txy_type[MAX_POINTS_PER_OBJECT];

        /* -----------------------------
        Variables
        --------------------------------*/

        database_trajectory DB_traj;
        database_nD DBn = new database_nD();

        int i, new_i, j, c, K_anonym, max_n, max_dim, n_equiv_classes, equiv_count;
        int i2, j2, k;
        int break_every_n, available_trash, trash_max, print_stats;
        double trash_max_percentage;
        double delta, semi_delta, delta_max, tmp_d, cluster_radius, largest_cluster_radius, max_translation;
        double[] object;
        double[] center;
        double TDC;
        int moved_points, moved_trajectories, one_point_moved;
        int trash_size, trash_size_as_points, too_small_classes;
        output_clusters_struct clusters;
        int DM_metric=0;
        double ID_metric;
        double CM_metric=0;

        /* Set mandatory parameters */
        K_anonym=anonymousK;
        delta=0;

        /* Default value for optional parameters */
        break_every_n=5;
        delta_max=0.01;
        trash_max_percentage=10.0;
        semi_delta=delta/2; /* Radius of the anonimity set tube */

        if(K_anonym<1) { printf("Error: parameter K must be > 0\n"); exit(-1); };
        if(delta<0) { printf("Error: parameter delta must be >= 0\n"); exit(-1); };
        if(break_every_n<1) { printf("Error: parameter break_every_n must be > 0\n"); exit(-1); };

        printf("Parameters:\n");
        printf(String.format("  K=%d, delta=%.3f, pi=%d, delta_max=%.3f, trash_max=%.1f%%\n",
                K_anonym, delta, break_every_n,delta_max, trash_max_percentage));

        /* -----------------------------
        Loading input data
        --------------------------------*/
        printf("Loading data...\n");

        DB_traj=load_objects_compact(txy_buffer, data, break_every_n);
        if(DB_traj.n==0) { printf("** The data file contains no objects! \n"); exit(-1); };
        /* Translate trash_max from percentage to an absolute value */
        trash_max=(int)(trash_max_percentage*(double)(DB_traj.n+DB_traj.traj_deleted)/100.0);

        /* -----------------------------
        Create equivalence classes
        --------------------------------*/
        printf("Creating equivalence classes...");
        /* : replace bubble sort with qsort or something else */ // Song:no use for we have equal length per traj
        Arrays.sort(DB_traj.db, new Comparator<trajectory>() {
            @Override
            public int compare(trajectory o1, trajectory o2) {
                if (o1.from > o2.from) {
                    return 1;
                } else {
                    if (o1.from == o2.from) {
                        return Integer.valueOf(o1.to).compareTo(o2.to);
                    } else {
                        return -1;
                    }
                }
            }
        });
//        for(i=DB_traj.n-1; i>0;)
//        {
//            new_i=0;
//            for(j=0; j<i; j++)
//                if( (DB_traj.db[j].from>DB_traj.db[j+1].from) ||
//                        ((DB_traj.db[j].from==DB_traj.db[j+1].from) && ((DB_traj.db[j].to>DB_traj.db[j+1].to))) )
//                {
//                    traj_tmp=DB_traj.db[j];
//                    DB_traj.db[j]=DB_traj.db[j+1];
//                    DB_traj.db[j+1]=traj_tmp;
//                    new_i=j;
//                }
//            i=new_i;
//        }
        printf("Done.\n");


        /* ---------------------------------------
        Process all equivalence classes
        ------------------------------------------*/
        printf("Cluster/translate/save...\r");

        /* Pre-scan to compute eq. classes & stats, and to make space for max-sized equivalence class */
        max_n=0;
        max_dim=0;
        n_equiv_classes=0;
        for(i=0; i<DB_traj.n; i=j)
        {
            for(j=i; (j<DB_traj.n) && (DB_traj.db[j].from==DB_traj.db[i].from) && (DB_traj.db[j].to==DB_traj.db[i].to); j++);
            max_n=Math.max(max_n, j-i);
            max_dim=Math.max(max_dim, 2*(DB_traj.db[i].to-DB_traj.db[i].from+1));
            n_equiv_classes++;
        }

        DBn.x= new double[max_n+1][]; /* "+1" for a dummy object --> trick in nwa.c */
        DBn.name= new int[max_n];
        center=new double[max_dim]; /* Dummy trajectory used in several place */

        /* Stats */
        TDC=0;
        largest_cluster_radius=-1;
        max_translation=-1;
        moved_points=0;
        moved_trajectories=0;
        trash_size=0;
        trash_size_as_points=0;
        available_trash=Math.max(0, trash_max-DB_traj.traj_deleted); /* Trajs. deleted in preprocessing count as trash */
        too_small_classes=0;

        /* Iterate for all equivalence classes */
        equiv_count=0;
        for(i=0; i<DB_traj.n; i=j)
        {
            equiv_count++;

            /* -----------------------------
            Clustering
            --------------------------------*/
            /* get first (i) and last (j) trajectory index of next equivalence class */
            for(j=i; (j<DB_traj.n) && (DB_traj.db[j].from==DB_traj.db[i].from) && (DB_traj.db[j].to==DB_traj.db[i].to); j++);

            /* Make DBn contain the data in DB_traj.db[i..j-1] */
            DBn.n=j-i;
            printf(String.format("Processing equivalence classes: %d of %d (%d trajs)           \r",
                    equiv_count, n_equiv_classes,DBn.n));
            DBn.dim=2*(DB_traj.db[i].to-DB_traj.db[i].from+1); /* =to-from+1 points */

            if(DBn.n >= K_anonym)  /* ...then we can anonymize this dataset */
            {
                for(k=i; k<j; k++)
                {
                    DBn.x[k-i]=DB_traj.db[k].x;
                    DBn.name[k-i]=DB_traj.db[k].name;
                }
                DBn.x[DBn.n]=center; /* Dummy element used as initial pivot (allocated once for all) */

                disFaster = new FasterEuclidean_nD(DBn);
                clusters = NwaParallel.never_walk_alone(DBn, K_anonym,
                        DB_traj.pseudo_diameter*delta_max,
                        (available_trash*DBn.n)/DB_traj.n); /* max trash is proportional to the size of eq. class */

                /* -----------------------------
                Space translation
                --------------------------------*/
                for(c=0; c<(clusters.n-1); c++)  /* scan all clusters (not trash) */
                {
                    /* Add contribution of the cluster to Discernibility Metric */
                    DM_metric += clusters.c_size[c] * clusters.c_size[c];   /* |cluster|^2 */

                    /* Compute the center of the cluster */
                    for(j2=0; j2<DBn.dim; j2++) center[j2]=0.0;
                    for(i2=0; i2<clusters.c_size[c]; i2++)
                        for(j2=0; j2<DBn.dim; j2++)
                            center[j2]+=DBn.x[ clusters.c_objects[c][i2] ][j2];
                    for(j2=0; j2<DBn.dim; j2++)
                        center[j2] = center[j2] / clusters.c_size[c];

                    /* move each object of the cluster towards the center (if needed) */
                    cluster_radius=0;
                    for(i2=0; i2<clusters.c_size[c]; i2++)
                    {
                        /* Update largest_cluster_radius = max distance between center and trajectories */
                        tmp_d=euclidean_nD(DBn, DBn.n, clusters.c_objects[c][i2]);
                        largest_cluster_radius=Math.max(largest_cluster_radius, tmp_d);
                        cluster_radius=Math.max(cluster_radius, tmp_d);
                        object=DBn.x[clusters.c_objects[c][i2]];
                        one_point_moved=0;
                        for(j2=0; j2<DBn.dim; j2+=2) /* Translation is done 2-by-2 coordinates, i.e., on 2D space */
                        {
                            tmp_d=Math.sqrt(SQUARE(object[j2]-center[j2])+SQUARE(object[j2+1]-center[j2+1]));
                            if(tmp_d>semi_delta) /* Translation is needed => do it directly in DBn */
                            {
                                /* Update max_translation */
                                max_translation=Math.max(max_translation, tmp_d-semi_delta);
                                /* Translate "object" towards "center" */
                                object[j2]=center[j2]+(object[j2]-center[j2])*semi_delta/tmp_d;
                                object[j2+1]=center[j2+1]+(object[j2+1]-center[j2+1])*semi_delta/tmp_d;
                                TDC+=tmp_d-semi_delta;
                                moved_points++;
                                one_point_moved=1;
                            }
                        }
                        if(one_point_moved==1) moved_trajectories++;
                    }
                    /* Update CM_metric */
                    CM_metric += (double)clusters.c_size[c]*cluster_radius;

                }
                trash_size+=clusters.c_size[clusters.n-1]; /* Last cluster contains the trash */
                trash_size_as_points+=clusters.c_size[clusters.n-1]*DBn.dim/2; /* Points trashed */

                /* -----------------------------
                Output D' (translated data)
                --------------------------------*/
//                for(c=0; c<(clusters.n-1);c++)  /* NOTE: last cluster (the trash) is not printed */ {
//                    for (i2 = 0; i2 < clusters.c_size[c]; i2++) {
//                        for (j2 = 0; j2 < DBn.dim; j2 += 2) {
//                            fprintf(OF, "%d\t%d\t%.6f\t%.6f\n",
//                                    DBn.name[clusters.c_objects[c][i2]],
//                                    (DB_traj.db[i].from + j2 / 2),
//                                    DBn.x[clusters.c_objects[c][i2]][j2],
//                                    DBn.x[clusters.c_objects[c][i2]][j2 + 1]);
//                        }
//                    }
//                }


                for(c=0; c<(clusters.n-1);c++)  /* NOTE: last cluster (the trash) is not printed */ {
                    for (i2 = 0; i2 < clusters.c_size[c]; i2++) {
                        List<GPSPoint> traj = new ArrayList<>();
                        for (j2 = 0; j2 < DBn.dim; j2 += 2) {
                            int id = DBn.name[clusters.c_objects[c][i2]];
                            int timeSlot = (DB_traj.db[i].from + j2 / 2);
                            double lon = DBn.x[clusters.c_objects[c][i2]][j2];
                            double lat = DBn.x[clusters.c_objects[c][i2]][j2 + 1];
                            GPSPoint p = new GPSPoint(id, timeSlot, lat, lon);
                            p.setTimeSlot(timeSlot);
                            traj.add(p);
                        }
                        result.add(traj);
                    }
                }


            } else /* ...we cannot cluster ==> do nothing, excepted updating statistics on trash*/
            {
                trash_size += DBn.n; /* Trajs. trashed */
                trash_size_as_points += DBn.n*DBn.dim/2; /* Points trashed */
                too_small_classes++;
            }

        }

        /* Finished ! */
        printf(String.format("Processing equivalence classes: Done!  [ %d equiv. classes ]               \n", n_equiv_classes));

        /* Optional statistics are printed on STDERR */
        /* Compute / finalize quality metrics */
        DM_metric += ((int)DB_traj.traj_deleted+(int)trash_size)*((int)DB_traj.n+(int)DB_traj.traj_deleted); /*   |trash|*|D|   */

        ID_metric = TDC + max_translation*(
                (double)DB_traj.points_deleted+(double)trash_size_as_points ); /*  TTD + Omega*pts_trash   */

        printf(String.format("%d %.5f %d %.5f %d ", K_anonym, delta, break_every_n, delta_max, trash_max)); /* parameters */
        printf(String.format("%d %d %d %d %.5f %d %d %d %d ", DB_traj.n, DB_traj.traj_deleted, DB_traj.n_points,
                DB_traj.points_deleted, DB_traj.pseudo_diameter,
                max_n, max_dim, n_equiv_classes, too_small_classes)); /* preprocessing */
        printf(String.format("%d %d %d %d %.5f %.5f %.5f %d %.5f %.5f\n", trash_size, trash_size_as_points, moved_trajectories,
                moved_points, TDC, largest_cluster_radius, max_translation, DM_metric, ID_metric, CM_metric)); /* clustering and translation */

        return result;

    }

    private static class NwaParallel {

        /**
         * NEVER_WALK_ALONE clustering algorithm
         */
        private static output_clusters_struct never_walk_alone(database_nD DB, int K, double first_max_radius, int first_max_trash_size) {
            /*  NOTICE: DB->x has to contain at least DB.n+1 vectors -- the extra one is a "dummy trajectory"
                The (DB.n+1)-th should be allocated as the others, though with any values */

            PriorityQueue<foint> myHeap; //smallest at top;
            foint myfoint;
            int i, S_size, pivot, best_pivot, j, last_cluster = 0, next_i;
            int[] assignments;
            int[] cluster_size;
            int[] alive;
            int[] pivots;
            byte[] ex_pivots;
            int n_ex_pivots;
            int[] k_buffer;
            double best_d, next_d, tmp_d, max_radius, farthest_d;
            int max_trash_size;
            double[] center;
            output_clusters_struct output_clusters = new output_clusters_struct();
            boolean solution_found;

            if (DB.n < K) {
                printf(String.format("ERROR: Never_walk_alone called with K (%d) greater than database_size (%d)!\n", K, DB.n));
                exit(-1);
            }

            /* Global Memory allocations and Initializations */

            /* Objects can be: (1) alive; (2) alive but already used as pivots - and failed; (3) dead */
            alive = new int[DB.n];       /* alive[0..S_size-1] contains IDs of alive objects */
            ex_pivots = new byte[DB.n]; /* ex_pivots[i]==0 <=> object "i" was never used as pivot */
            pivots = new int[DB.n];      /* pivots[i] = pivot/center of cluster "i" */
            k_buffer = new int[K];
            assignments = new int[DB.n]; /* assignments[i] = cluster of object "i" */
            cluster_size = new int[DB.n]; /* (oversized, to avoid problems...) */

            myHeap = new PriorityQueue<>();//HeapAlloc(K,heap_cmp_reverse); /* Heap used for efficiently collecting (K-1)-closest objects */

            /* First pivot = center of the dataset --> trick: uses a dummy object of the dataset */
            center = DB.x[DB.n];
            for (j = 0; j < DB.dim; j++) {
                center[j] = 0.0;
            }
            for (i = 0; i < DB.n; i++) {
                for (j = 0; j < DB.dim; j++) {
                    center[j] += DB.x[i][j];
                }
            }
            for (j = 0; j < DB.dim; j++) {
                center[j] = center[j] / DB.n;
            }

            /* Cluster constraints - starting values*/
            max_radius = first_max_radius * Math.sqrt((double) DB.dim / 2.0); /* Takes into account dimensionality */
            max_trash_size = first_max_trash_size;

            /* Here the external cycle starts */
            solution_found = false;
            while (!solution_found) {
                /* printf(" [max_radius=%f] ",max_radius); fflush(NULL); */
                /* Init variables */
                for (i = 0; i < DB.n; i++) {
                    alive[i] = i;
                    ex_pivots[i] = 0;
                    cluster_size[i] = 0; /* CHECK: not needed */
                }
                S_size = DB.n;   /* Number of alive objects */
                n_ex_pivots = 0;  /* Number of alive yet ex-pivot objects */
                last_cluster = -1;
                next_i = -1;	/* Farthest pivot from actual one (-1 = unknown) */
                next_d = -1;
                pivot = DB.n;	/* First pivot = center of the dataset */

                /* Step 1: find clusters of size K and radius <= max_radius */
                while ((S_size >= K) && (S_size > n_ex_pivots)) { /* Stop if not enough objects alive or all are ex-pivots */

                    last_cluster++;
                    cluster_size[last_cluster] = 0;

                    /* Choose next pivot -- excluding ex-pivots, that cannot be pivot again */
                    if (next_i != -1) {
                        i = next_i;
                        best_d = next_d;
                    } else {
                        best_d = -1;
                        i = -1;
                        for (j = 0; j < S_size; j++) /* scan all alive objects */
                            if (ex_pivots[alive[j]] == 0 && ((tmp_d = euclidean_nD(DB, pivot, alive[j])) > best_d)) {
                                i = j;
                                best_d = tmp_d;
                            }
                    }
                    if (i == -1) {
                        printf("ERROR: inconsistency looking for next pivot!\n");
                        exit(-1);
                    }
                    pivot = alive[i]; /* Careful: "i" has to be keeped till the end of the cycle! */

                    /* Select closest K-1 objects and find next possible pivot */
                    next_d = -1;
                    next_i = -1;
                    for (j = 0; j < S_size; j++) {
                        if (alive[j] != pivot) {
                            myHeap.add(new foint(j, euclidean_nD(DB, pivot, alive[j])));
                            if (myHeap.size() >= K) /* keep only top K-1 in the heap */ {
                                myfoint = myHeap.poll();
                                /* Check whether the object can be next pivot */
                                if ((myfoint.v.d > next_d) && (ex_pivots[alive[myfoint.v.i]] == 0)) {
                                    next_i = myfoint.v.i;
                                    next_d = myfoint.v.d;
                                }
                                //                            free(myfoint.v);    /* Remove farthest object */
                            }
                        }
                    }

                    /* Get the distance of farthest object -- it's on the top of the heap */
                    if (myHeap.size() > 0) {
                        myfoint = myHeap.peek();
                        farthest_d = myfoint.v.d;
                    } else
                        farthest_d = 0; /* Very strange case: occurs only when S_size = K = 1 */

                    /* Pop all the content of the heap, and copy the "i"s (indexes of alive[]) to k_buffer */
                    for (j = 0; j < K - 1; j++) {
                        myfoint = myHeap.poll();
                        k_buffer[j] = myfoint.v.i;
                    }

                    /* Check cluster diameter (more exactly: max distance from cluster pivot) */
                    if (farthest_d <= max_radius) /* Check cluster radius (="pseudo-radius") */ {
                        /* 1. "pivot" becomes a real pivot */
                        pivots[last_cluster] = pivot;

                        /* 2. From now on, check out all the members of the cluster */
                        /* Add pivot to k_buffer, whose elements are going to be "checked out" */
                        k_buffer[K - 1] = i; /* "i" is the position of the pivot we are deleting */

                        /* To avoid problems with alive[], remove objects from right to left */
                        /* OPTIMIZE: can we avoid sorting? */
                        Arrays.sort(k_buffer);
                        for (j = K - 1; j >= 0; j--) {
                            assignments[alive[k_buffer[j]]] = last_cluster;
                            cluster_size[last_cluster]++;
                            if (ex_pivots[alive[k_buffer[j]]] == 1) /* remove dead ex-pivots */
                                n_ex_pivots--;
                            /* Virtually pop out the object from alive[] */
                            alive[k_buffer[j]] = alive[S_size - 1];
                            S_size--;
                            if (next_i == S_size) /* I moved my next pivot... */
                                next_i = k_buffer[j];
                        }
                    } else {
                        /* do nothing, just "deactivate" the pivot */
                        ex_pivots[pivot] = 1;
                        n_ex_pivots++;
                        last_cluster--; /* undo last step */
                    }
                }

                /* Step 2: assign unassigned objects to existing clusters (where possible) */
                /*         alive[] contains the remaining objects: try to assign them  */
                for (i = 0; i < S_size; i++) {
                    if (last_cluster >= 0) /* If there is at least one cluster, find closest one */
                    {
                        best_pivot = 0;
                        best_d = euclidean_nD(DB, pivots[0], alive[i]);
                        for (j = 1; j <= last_cluster; j++) {
                            if ((tmp_d = euclidean_nD(DB, pivots[j], alive[i])) < best_d) {
                                best_pivot = j;
                                best_d = tmp_d;
                            }
                        }
                        if (best_d > max_radius)
                            best_pivot = last_cluster + 1; /* This is the trash */
                    } else
                        best_pivot = last_cluster + 1; /* This is the trash */

                    assignments[alive[i]] = best_pivot;
                    cluster_size[best_pivot]++;
                }

                /* Check trash size */
                solution_found = (cluster_size[last_cluster + 1] <= max_trash_size);
                if (!solution_found) {
                    max_radius = max_radius * 1.5;		/* Exponential-growth radius threshold */
                    max_trash_size = max_trash_size;	/* Constant-valued trash threshold */
                    printf(String.format("solution not found, increase max_radius to %f", max_radius));
                }
            }

            /* Package results for output */
            output_clusters.n = last_cluster + 2;
            output_clusters.c_size = cluster_size;
            output_clusters.c_objects = new int[output_clusters.n][];
            for (i = 0; i < output_clusters.n; i++)
                output_clusters.c_objects[i] = new int[cluster_size[i]];
            for (i = 0; i < output_clusters.n; i++)
                cluster_size[i] = 0;
            for (i = 0; i < DB.n; i++) {
                output_clusters.c_objects[assignments[i]][cluster_size[assignments[i]]] = i;
                cluster_size[assignments[i]]++;
            }
            return output_clusters;
        }
    }

    private static class FasterEuclidean_nD
    {
        private double[][] memory;

        FasterEuclidean_nD(final database_nD db) {
            memory = new double[db.n][];
            for(int i=0; i<db.n; i++){
                memory[i] = new double[db.n];
            }
            for(int i=0; i<db.n-1; i++)
            {
                for (int j = i + 1; j < db.n; j++) {
                    double dis = original_euclidean_nD(db, i, j);
                    memory[i][j] = dis;
                    memory[j][i] = dis;
                }
            }
        }

        void FasterEuclidean_nD_Parallel(final database_nD db) {
            memory = new double[db.n][];
            for(int i=0; i<db.n; i++){
                memory[i] = new double[db.n];
            }
            int parallelCount = 2;
            final Semaphore semaphore = new Semaphore(parallelCount);
            List<Thread> threadList = new ArrayList<>(db.n-1);
            for(int i=0; i<db.n-1; i++)
            {
                final int finalI = i;
                Thread t = new Thread(new Runnable() {
                    @Override
                    public void run() {
                        try {
                            semaphore.acquire();
                            for (int j = finalI + 1; j < db.n; j++) {
                                double dis = original_euclidean_nD(db, finalI, j);
                                memory[finalI][j] = dis;
                                memory[j][finalI] = dis;
                            }
                        } catch (InterruptedException e) {
                            e.printStackTrace();
                        } finally {
                            semaphore.release();
                        }
                    }
                });
                threadList.add(t);
                t.start();
            }
            try {
                for(Thread t: threadList) {
                    t.join();
                }
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
        }

        public Double get(int a, int b)
        {
            return memory[a][b];
        }

        public void verify(int a, int b, double dis) {
            if(memory[a][b]!=dis){
                throw new RuntimeException("dis not equal! "+dis+" "+memory[a][b]);
            }
        }
    }


    public static class GPSPoint extends GPXEntry
    {
        private int carID;

        public GPSPoint(int carID, long timeSlot, double longitude, double latitude)
        {
            super(latitude, longitude, timeSlot);
            this.carID = carID;
        }

        public int getCarID()
        {
            return carID;
        }

        public void setCarID(int carID)
        {
            this.carID = carID;
        }

        public int getTimeSlot()
        {
            if(this.getTime()<Integer.MAX_VALUE) return (int) this.getTime();
            else throw new RuntimeException("timestamp larger than INT.max_value!");
        }

        public void setTimeSlot(int timeSlot)
        {
            this.setTime(timeSlot);
        }

        public String toStringShort(char sep){
            return String.valueOf(carID) + sep +
                    getTime() + sep +
                    getLat() + sep +
                    getLon();
        }

        public String toString(){
            return toStringShort(' ');
        }

    }

}
