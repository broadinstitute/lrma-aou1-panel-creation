import java.util.Arrays;
import java.util.Random;
import java.util.zip.GZIPInputStream;

import java.io.IOException;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.BufferedWriter;

import java.awt.Color;
import java.awt.Stroke;
import java.awt.BasicStroke;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;


/**
 * Assume that we are given a multi-sample VCF with SVs and SNPs. Consider a 
 * maximal set of calls W that might overlap in the reference. The genotypes of
 * such calls in a sample S might be such that some calls overlap on the same 
 * haplotype of S. The program finds a max-weight set of non-overlapping 
 * calls for each W, S and haplotype, and it modifies the GTs accordingly.
 *
 * This is useful if downstream scripts remove all the calls in W that occur on 
 * a haplotype of S, if there is even a single pair of overlapping calls on such
 * a haplotype. See e.g.:
 * 
 * https://bitbucket.org/jana_ebler/vcf-merging/src/cb56c5df703cad11e7b777623d15a0b839931685/pangenome-graph-from-callset/scripts/merge_vcfs.py#lines-419
 */
public class FixVariantCollisions {
    /**
     * Weight tag (in the INFO field or in the SAMPLE field of the VCF).
     */
    private static String WEIGHT_TAG;
    private static boolean WEIGHT_TAG_IN_SAMPLE_COLUMN;
    
    /**
     * FALSE=remove entire VCF records from a sample;
     * TRUE=remove single ones from a sample's GT.
     */
    private static boolean METHOD;
    
    /**
     * Variant types
     */
    private static final int DEL = 0;
    private static final int INV = 1;
    private static final int DUP = 2;
    private static final int INS = 3;
    private static final int SNP = 4;
    private static final int REPLACEMENT = 5;
    
    /**
     * Genotypes
     */
    private static final int PHASED_00 = 0;
    private static final int PHASED_01 = 1;
    private static final int PHASED_10 = 2;
    private static final int PHASED_11 = 3;
    private static final int UNPHASED_00 = 4;
    private static final int UNPHASED_01 = 5;
    private static final int UNPHASED_10 = 6;
    private static final int UNPHASED_11 = 7;
    
    private static final int PHASED_D0 = 8;
    private static final int PHASED_0D = 9;
    private static final int PHASED_D1 = 10;
    private static final int PHASED_1D = 11;
    private static final int PHASED_DD = 12;
    private static final int UNPHASED_D0 = 13;
    private static final int UNPHASED_0D = 14;
    private static final int UNPHASED_D1 = 15;
    private static final int UNPHASED_1D = 16;
    private static final int UNPHASED_DD = 17;
    
    /**
     * A maximal set of overlapping intervals.
     *
     * Remark: the last interval in $window$ might not overlap with the previous
     * intervals, but it might be the first interval of the next set of 
     * overlapping intervals.
     */
    private static Interval[] window;
    
    /**
     * Last element in $window$.
     */
    private static int windowLast;
    
    /**
     * Last position of an overlapping interval in $window$ (one-based,
     * inclusive). Might be smaller than the last position of the last interval
     * in $window$, if such an interval does not overlap with the previous ones.
     */
    private static int windowLastPos;
    

    /**
     * Example usage:
     *
     * java FixVariantCollisions input.vcf.gz 0 WEIGHT_TAG 0 output.vcf windows.txt histogram.txt null
     *
     * @param args
     * Input arguments:
     * 0: input VCF.GZ; every record is assumed to be biallelic;
     * 1: operations allowed to fix the genotypes of a sample: 0=can only remove
     *    an entire VCF record; 1=can remove single ones from a GT;
     * 2: ID of the weight field; if this field is not found, all weights are
     *    set to one; weights are assumed to be non-negative;
     * 3: given a VCF record in a sample, assign it a weight encoded in the
     *    sample column (1) or in the INFO field (0).
     *
     * Output arguments:
     * 4: writes to this file (uncompressed VCF) the input VCF.GZ with all
     *    collisions fixed; "null"=collisions are not fixed and the file is not
     *    written;
     * 5: writes to this file the list of all windows processed;
     * 6: histogram: for each X, the number of (window,sample) pairs with X
     *    collisions; "null"=the histogram is not printed;
     * 7: directory where to store figures for every window and sample with
     *    collisions; "null"=figures are not printed. 
     */
    public static void main(String[] args) throws IOException {
        // Input arguments
        final String INPUT_VCF_GZ = args[0];
        METHOD=args[1].equalsIgnoreCase("1");
        WEIGHT_TAG=args[2];
        WEIGHT_TAG_IN_SAMPLE_COLUMN=Integer.parseInt(args[3])==1;
        // Output arguments
        final String OUTPUT_VCF = args[4];
        final String OUTPUT_WINDOWS = args[5];
        final String OUTPUT_HISTOGRAM = args[6];
        String OUTPUT_FIGURES_DIR = args[7];
        
        int i;
        int nRecords;
        BufferedReader brVCF;
        BufferedWriter bwVCF, bwWindows, bwHistogram;
        long[] histogram;
        
        initType2color();
        window = new Interval[100];
        for (i=0; i<window.length; i++) window[i] = new Interval();
        histogram = new long[100];
        windowLast=-1; windowLastPos=-1;
        brVCF = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(INPUT_VCF_GZ))));
        if (OUTPUT_VCF.equalsIgnoreCase("null")) bwVCF=null;
        else bwVCF = new BufferedWriter(new FileWriter(OUTPUT_VCF));
        if (OUTPUT_WINDOWS.equalsIgnoreCase("null")) bwWindows=null;
        else bwWindows = new BufferedWriter(new FileWriter(OUTPUT_WINDOWS));
        if (OUTPUT_FIGURES_DIR.equalsIgnoreCase("null")) OUTPUT_FIGURES_DIR=null;
        nRecords=0;
        while (true) {
            loadNextWindow(brVCF,bwVCF);
            if (windowLast==-1) break;
            fixWindow(bwVCF,bwWindows,OUTPUT_FIGURES_DIR,histogram);
            nRecords+=(isLastInWindow()?windowLast:windowLast-1)+1;
            if (nRecords%10000==0) System.err.println("Loaded "+nRecords+" records");
        }
        brVCF.close(); bwVCF.close(); bwWindows.close();
        if (!OUTPUT_HISTOGRAM.equalsIgnoreCase("null")) {
            bwHistogram = new BufferedWriter(new FileWriter(OUTPUT_HISTOGRAM));
            bwHistogram.write("#nCollision \t nHaplotypes\n");
            for (i=0; i<histogram.length; i++) bwHistogram.write(i+"\t"+histogram[i]+"\n");
            bwHistogram.close();
        }
    }
    
    
    /**
     * Loads in $window$ a maximal set of overlapping intervals, regardless of
     * their genotypes.
     *
     * Remark: after this procedure completes, the last interval in $window$
     * might not overlap with the previous ones.
     *
     * @param bw if not null, input headers (if any) are written here.
     */
    private static final void loadNextWindow(BufferedReader br, BufferedWriter bw) throws IOException {
        int i;
        String str;
        Interval tmpInterval;
        
        if (isLastInWindow()) {
            windowLast=-1;
            windowLastPos=-1;
        }
        else {
            tmpInterval=window[0];
            window[0]=window[windowLast];
            window[windowLast]=tmpInterval;
            windowLastPos=window[0].last;
            windowLast=0;
        }
        while (true) {
            str=br.readLine();
            if (str==null) break;
            if (str.charAt(0)==COMMENT) {
                if (bw!=null) { bw.write(str); bw.newLine(); }
                continue;
            }
            windowLast++;
            if (windowLast==window.length) {
                Interval[] newArray = new Interval[window.length<<1];
                System.arraycopy(window,0,newArray,0,window.length);
                for (i=window.length; i<newArray.length; i++) newArray[i] = new Interval();
                window=newArray;
            }
            tmpInterval=window[windowLast];
            tmpInterval.init(str);
            if (!isLastInWindow()) break;
            if (tmpInterval.last>windowLastPos) windowLastPos=tmpInterval.last;
        }
        for (i=0; i<=windowLast; i++) window[i].inputIndex=i;
    }
    
    
    /**
     * @return TRUE iff the last interval of $window$ overlaps with the previous
     * intervals.
     */
    private static final boolean isLastInWindow() {
        return windowLast==-1 || windowLastPos==-1 || (window[windowLast].chr==window[0].chr && window[windowLast].first<=windowLastPos);
    }
    
    
    /**
     * @param bwVCF if NULL, the procedure does not fix collisions;
     * @param figuresDir if not NULL, the procedure stores in this directory an
     * overlap diagram of the window for every sample that contains a collision
     * (one PNG file per sample);
     * @param histogram for each $i$, the number of (window,sample) pairs with
     * $i$ collisions; windows with exactly one VCF record are not considered in
     * the count.
     */
    private static final void fixWindow(BufferedWriter bwVCF, BufferedWriter bwWindows, String figuresDir, long[] histogram) throws IOException {
        final int FIRST = window[0].first;
        final int N_COLUMNS = (2+windowLastPos-FIRST+1)*PIXELS_PER_POS;
        final int LAST = isLastInWindow()?windowLast:windowLast-1;
        final int N_SAMPLES = window[0].genotypes.length;
        final int CHR = window[0].chr;
        final String WINDOW_DIR = figuresDir==null?null:figuresDir+"/chr"+CHR+"_"+FIRST+"_"+windowLastPos;
        
        int i, j;
        int nCollisions;
        Random random = null;
        File directory = null;
        BufferedImage image = null;
        Graphics2D graphics = null;
        
        if (bwWindows!=null) bwWindows.write(window[0].chr+"\t"+window[0].first+"\t"+windowLastPos+"\t"+(LAST+1));
        if (LAST==0) { 
            if (bwWindows!=null) {
                for (i=0; i<N_SAMPLES; i++) bwWindows.write("\t0");
                bwWindows.newLine();
            }
            if (bwVCF!=null) window[0].toVCF(bwVCF); 
            return; 
        }
        Interval.order=Interval.ORDER_LAST_POS;
        Arrays.sort(window,0,LAST+1);
        if (figuresDir!=null) {
            random = new Random();
            directory = new File(WINDOW_DIR);
            image = new BufferedImage(N_COLUMNS,N_ROWS,BufferedImage.TYPE_INT_RGB);
            graphics=image.createGraphics();
        }
        for (j=0; j<N_SAMPLES; j++) {
            nCollisions=countCollisions(LAST,j);
            histogram[nCollisions>histogram.length-1?histogram.length-1:nCollisions]++;
            if (bwWindows!=null) bwWindows.write("\t"+nCollisions);
            if (nCollisions==0) continue;
            if (figuresDir!=null) {
                if (!directory.exists()) directory.mkdirs();
                drawWindow(j,LAST,graphics,N_COLUMNS,FIRST,random);
                ImageIO.write(image,"png",new File(WINDOW_DIR+"/sample"+j+"_before.png"));
            }
            if (bwVCF!=null) {
                if (METHOD) independentSet2(j,LAST);
                else independentSet1(j,LAST);
                if (figuresDir!=null) {
                    drawWindow(j,LAST,graphics,N_COLUMNS,FIRST,random);
                    ImageIO.write(image,"png",new File(WINDOW_DIR+"/sample"+j+"_after.png"));
                }
            }
        }
        bwWindows.newLine();
        if (bwVCF!=null) {
            Interval.order=Interval.ORDER_INPUT;
            Arrays.sort(window,0,LAST+1);
            for (i=0; i<=LAST; i++) window[i].toVCF(bwVCF);
        }
    }
    
    
    /**
     * Remark: the procedure assumes that $window$ is sorted by last position.
     *
     * @return the number of pairs of intervals in $window[0..last]$ that 
     * collide in $sample$, assuming that a diploid call consists of two 
     * intervals.
     */
    private static final int countCollisions(int last, int sample) {
        int i, j;
        int startI, endI, gtI, out;
        
        out=0;
        for (i=last; i>=0; i--) {
            if (!window[i].isPresent(sample)) continue;
            startI=window[i].first; gtI=window[i].genotypes[sample];
            for (j=i-1; j>=0; j--) {
                if (window[j].last<startI) break;
                out+=N_GT_COLLISIONS[window[j].genotypes[sample]][gtI];
            }
        }
        return out;
    }
    
    
    
    
    // ------------------------ GENOTYPE PROCEDURES ----------------------------
    
    /**
     * @param gt a biallelic genotype;
     * @return a unique ID for every possible biallelic genotype.
     */
    private static final int gt2Id(String gt) {
        final char a = gt.charAt(0);
        final char b = gt.charAt(1);
        final char c = gt.charAt(2);
    
        if (b=='/') {
            if (a=='.') {
                if (c=='.') return UNPHASED_DD;
                else if (c=='0') return UNPHASED_D0;
                else return UNPHASED_D1;
            }
            else if (a=='0') {
                if (c=='.') return UNPHASED_0D;
                else if (c=='0') return UNPHASED_00;
                else return UNPHASED_01;
            }
            else {
                if (c=='.') return UNPHASED_1D;
                else if (c=='0') return UNPHASED_10;
                else return UNPHASED_11;
            }
        }
        else {
            if (a=='.') {
                if (c=='.') return PHASED_DD;
                else if (c=='0') return PHASED_D0;
                else return PHASED_D1;
            }
            else if (a=='0') {
                if (c=='.') return PHASED_0D;
                else if (c=='0') return PHASED_00;
                else return PHASED_01;
            }
            else {
                if (c=='.') return PHASED_1D;
                else if (c=='0') return PHASED_10;
                else return PHASED_11;
            }
        }
    }
    
    
    private static final String[] GT2STR = new String[] {
        "0|0", "0|1", "1|0", "1|1",  
        "0/0", "0/1", "1/0", "1/1",
        ".|0", "0|.", ".|1", "1|.", ".|.",  
        "./0", "0/.", "./1", "1/.", "./."
    };
    
    
    /**
     * Rows, columns: biallelic GT ids. Cells: {0,1,2} the number of pairs of
     * intervals that collide, assuming that a diploid call consists of two
     * intervals.
     */
    private static final int[][] N_GT_COLLISIONS = new int[][] { 
        // 0|0, 0|1, 1|0, 1|1,  
        // 0/0, 0/1, 1/0, 1/1,  
        // .|0, 0|., .|1, 1|., .|.,  
        // ./0, 0/., ./1, 1/., ./.
        {0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0},
        {0,1,0,1, 0,0,0,1, 0,0,1,0,0, 0,0,0,0,0},
        {0,0,1,1, 0,0,0,1, 0,0,0,1,0, 0,0,0,0,0},
        {0,1,1,2, 0,0,0,2, 0,0,1,1,0, 0,0,1,1,0},
        {0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0},
        {0,0,0,1, 0,0,0,1, 0,0,0,0,0, 0,0,0,0,0},
        {0,0,0,1, 0,0,0,1, 0,0,0,0,0, 0,0,0,0,0},
        {0,1,1,2, 0,0,0,2, 0,0,1,1,0, 0,0,1,1,0},
            
        {0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0},
        {0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0},
        {0,1,0,1, 0,0,0,1, 0,0,1,0,0, 0,0,0,0,0},
        {0,0,1,1, 0,0,0,1, 0,0,0,1,0, 0,0,0,0,0},
        {0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0},
        {0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0},
        {0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0},
        {0,0,0,1, 0,0,0,1, 0,0,0,0,0, 0,0,0,0,0},
        {0,0,0,1, 0,0,0,1, 0,0,0,0,0, 0,0,0,0,0},
        {0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0}
    };

    
    /**
     * Replaces a one with a zero on hap1.
     */
    private static final int[] REMOVE_HAP1_0 = new int[] {
        // 0|0, 0|1, 1|0, 1|1,
        // 0/0, 0/1, 1/0, 1/1,
        // .|0, 0|., .|1, 1|., .|.,
        // ./0, 0/., ./1, 1/., ./.
        PHASED_00,PHASED_01,PHASED_00,PHASED_01,  
        UNPHASED_00,UNPHASED_01,UNPHASED_10,UNPHASED_01,
        PHASED_D0,PHASED_0D,PHASED_D1,PHASED_0D,PHASED_DD,
        UNPHASED_D0,UNPHASED_0D,UNPHASED_D1,UNPHASED_1D,UNPHASED_DD
    };
    
    
    /**
     * Replaces a one with a dot on hap1.
     */
    private static final int[] REMOVE_HAP1_D = new int[] {
        // 0|0, 0|1, 1|0, 1|1,
        // 0/0, 0/1, 1/0, 1/1,
        // .|0, 0|., .|1, 1|., .|.,
        // ./0, 0/., ./1, 1/., ./.
        PHASED_00,PHASED_01,PHASED_D0,PHASED_D1,
        UNPHASED_00,UNPHASED_01,UNPHASED_10,UNPHASED_D1,
        PHASED_D0,PHASED_0D,PHASED_D1,PHASED_DD,PHASED_DD,
        UNPHASED_D0,UNPHASED_0D,UNPHASED_D1,UNPHASED_1D,UNPHASED_DD
    };
    
    
    /**
     * Replaces a one with a zero on hap2.
     */
    private static final int[] REMOVE_HAP2_0 = new int[] {
        // 0|0, 0|1, 1|0, 1|1,
        // 0/0, 0/1, 1/0, 1/1,
        // .|0, 0|., .|1, 1|., .|.,
        // ./0, 0/., ./1, 1/., ./.
        PHASED_00,PHASED_00,PHASED_10,PHASED_10,
        UNPHASED_00,UNPHASED_01,UNPHASED_10,UNPHASED_01,
        PHASED_D0,PHASED_0D,PHASED_D0,PHASED_1D,PHASED_DD,
        UNPHASED_D0,UNPHASED_0D,UNPHASED_D1,UNPHASED_1D,UNPHASED_DD
    };
    
    
    /**
     * Replaces a one with a dot on hap2.
     */
    private static final int[] REMOVE_HAP2_D = new int[] {
        // 0|0, 0|1, 1|0, 1|1,
        // 0/0, 0/1, 1/0, 1/1,
        // .|0, 0|., .|1, 1|., .|.,
        // ./0, 0/., ./1, 1/., ./.
        PHASED_00,PHASED_0D,PHASED_10,PHASED_1D,
        UNPHASED_00,UNPHASED_01,UNPHASED_10,UNPHASED_D1,
        PHASED_D0,PHASED_0D,PHASED_DD,PHASED_1D,PHASED_DD,
        UNPHASED_D0,UNPHASED_0D,UNPHASED_D1,UNPHASED_1D,UNPHASED_DD
    };
    
    
    

    // --------------------- INDEPENDENT SET PROCEDURES ------------------------
    
    /**
     * The optimal solution
     */
    private static double maxWeight;
    
    /**
     * Reused space
     */
    private static Interval[] sampleWindow;
    private static int sampleWindowLast;
    private static Point[] endpoints;
    private static int lastEndpoint;
    
    
    /**
     * Assume that we can fix the genotypes of $sample$ only by deleting records
     * on both haplotypes. Then we would like to keep a heaviest set of VCF 
     * records that do not overlap on any haplotype. The procedure computes a 
     * max-weight independent set of the graph where every element of 
     * $window[0..last]$ that occurs in $sample$ (on any haplotype) is a node, 
     * and where there is an edge iff two elements overlap on some haplotype.
     *
     * Remark: the procedure considers also unphased calls, since they might
     * collide with phased calls.
     *
     * Remark: the algorithm is a simple tree search on the (exponential) space
     * of all possible solutions, with a simple pruning based on the current
     * best solution.
     * 
     * Remark: the procedure assumes that $window[0..last]$ is sorted by last
     * position.
     */
    private static final void independentSet1(int sample, int last) {
        boolean error;
        int i, j;
        double weight, currentWeight;
        Interval currentInterval, neighbor;
        boolean[] active;

        // Collecting only the intervals that are present in $sample$ (some of
        // which might be unphased), and assigning a weight to them.
        if (sampleWindow==null || sampleWindow.length<last+1) sampleWindow = new Interval[last+1];
        sampleWindowLast=-1;
        for (i=0; i<=last; i++) {
            currentInterval=window[i];
            if (!currentInterval.isPresent(sample)) continue;
            sampleWindow[++sampleWindowLast]=currentInterval;
            currentInterval.clearIndependentSetVariables();
            currentInterval.setWeight(sample);
        }
        
        // Recursion
        maxWeight=0;
        active = new boolean[sampleWindowLast+1];
        Arrays.fill(active,true);
        for (i=sampleWindowLast; i>=0; i--) independendSet1_impl(i,active,new int[0],0,sample);
        error=false;
        if (maxWeight==0) error=true;
        else {
            error=true;
            for (i=0; i<=sampleWindowLast; i++) {
                if (sampleWindow[i].inIndependentSet) { error=false; break; }
            }
        }
        if (error) {
            System.err.println("No selected interval in the optimal solution of sample "+sample+"?!  maxWeight="+maxWeight+" sampleWindowLast="+sampleWindowLast);        
            System.exit(1);
        }
        
        // Removing from $sample$ every interval that is not in the max-weight
        // independent set.
        markIndependentSetOverlaps(sample,true,true,sampleWindowLast);
        for (i=0; i<=sampleWindowLast; i++) {
            if (sampleWindow[i].inIndependentSet) continue;
            if (sampleWindow[i].overlapsIS_hap1) sampleWindow[i].genotypes[sample]=REMOVE_HAP1_D[sampleWindow[i].genotypes[sample]];
            else sampleWindow[i].genotypes[sample]=REMOVE_HAP1_0[sampleWindow[i].genotypes[sample]];
            if (sampleWindow[i].overlapsIS_hap2) sampleWindow[i].genotypes[sample]=REMOVE_HAP2_D[sampleWindow[i].genotypes[sample]];
            else sampleWindow[i].genotypes[sample]=REMOVE_HAP2_0[sampleWindow[i].genotypes[sample]];
        }
    }
    
    
    /**
     * @param id position in $sampleWindow$ to solve for;
     * @param active one flag per position in $sampleWindow$; indicates whether 
     * the corresponding interval is active before solving for $id$;
     * @param ancestors ancestors of $id$ in the search tree;
     * @param ancestorsWeight sum of weights of all the elements in $ancestors$.
     */
    private static final void independendSet1_impl(int id, boolean[] active, int[] ancestors, double ancestorsWeight, int sample) {
        final Interval currentInterval = sampleWindow[id];
        final double weight_prime = ancestorsWeight+currentInterval.weight;
        
        boolean hasChild;
        int i;
        double upperBound;
        boolean[] active_prime;
        int[] ancestors_prime;
        
        active_prime = new boolean[active.length];
        System.arraycopy(active,0,active_prime,0,active.length);
        active_prime[id]=false;
        hasChild=false; upperBound=weight_prime;
        for (i=id-1; i>=0; i--) {
            if (active_prime[i] && !sampleWindow[i].precedes(currentInterval,sample)) active_prime[i]=false;
            if (active_prime[i]) {
                hasChild=true;
                upperBound+=sampleWindow[i].weight;
            }
        }
        if (upperBound<maxWeight) return;  // Pruning all recursive calls
        if (hasChild) {
            ancestors_prime = new int[ancestors.length+1];
            System.arraycopy(ancestors,0,ancestors_prime,0,ancestors.length);
            ancestors_prime[ancestors_prime.length-1]=id;
            for (i=id-1; i>=0; i--) {
                if (!active_prime[i]) continue;
                if (upperBound<maxWeight) break;  // Pruning this and all the following recursive calls
                independendSet1_impl(i,active_prime,ancestors_prime,weight_prime,sample);
                upperBound-=sampleWindow[i].weight;
            }
        }
        else if (weight_prime>maxWeight) {  // Base case of the recursion
            maxWeight=weight_prime;
            for (i=0; i<=sampleWindowLast; i++) sampleWindow[i].inIndependentSet=false;
            for (i=0; i<ancestors.length; i++) sampleWindow[ancestors[i]].inIndependentSet=true;
            currentInterval.inIndependentSet=true;
        }
    }
    
    
    /**
     * Assume that we can fix the genotypes of $sample$ only by deleting
     * specific calls on specific haplotypes. Then we would like to keep, on
     * each haplotype, a heaviest set of calls that do not overlap on that
     * haplotype. The procedure computes, for each haplotype, a max-weight
     * independent set of the graph where every element of $window[0..last]$
     * that occurs on that haplotype is a node, and where there is an edge iff
     * two elements overlap.
     *
     * Remark: the procedure considers only phased calls and 1/1 unphased calls.
     *
     * Remark: the procedure assumes that $window[0..last]$ is sorted by last
     * position.
     */
    private static final void independentSet2(int sample, int last) {
        int i;
        
        // Allocating reused space
        if (sampleWindow==null || sampleWindow.length<last+1) sampleWindow = new Interval[last+1];
        if (endpoints==null) {
            endpoints = new Point[(last+1)<<1];
            for (i=0; i<endpoints.length; i++) endpoints[i] = new Point();
        }
        else if (endpoints.length<(last+1)<<1) {
            Point[] newArray = new Point[(last+1)<<1];
            System.arraycopy(endpoints,0,newArray,0,endpoints.length);
            for (i=endpoints.length; i<newArray.length; i++) newArray[i] = new Point();
            endpoints=newArray;
        }
        
        // Solving Hap 1
        sampleWindowLast=-1;
        for (i=0; i<=last; i++) {
            if (window[i].onHap1(sample)) {
                sampleWindow[++sampleWindowLast]=window[i];
                window[i].setWeight(sample);
            }
        }
        if (sampleWindowLast!=-1) independentSet2_impl(sample,true);
        
        // Solving Hap 2
        sampleWindowLast=-1;
        for (i=0; i<=last; i++) {
            if (window[i].onHap2(sample)) {
                sampleWindow[++sampleWindowLast]=window[i];
                window[i].setWeight(sample);
            }
        }
        if (sampleWindowLast!=-1) independentSet2_impl(sample,false);
    }


    /**
     * Implements the linear algorithm in Section 2 of:
     *
     * Hsiao et al. "An efficient algorithm for finding a maximum weight
     * 2-independent set on interval graphs." Information Processing Letters
     * 43.5 (1992): 229-235.
     *
     * @param onHap1 the procedure assumes that all intervals in $sampleWindow$ 
     * occur on haplotype 1 (TRUE) or on haplotype 2 (FALSE) of $sample$.
     */
    private static final void independentSet2_impl(int sample, boolean onHap1) {
        boolean error;
        int i, j;
        double availableWeight;
        Interval currentInterval, previous;

        // Computing a max-weight independent set
        intervals2endpoints();
        for (i=0; i<=sampleWindowLast; i++) sampleWindow[i].clearIndependentSetVariables();
        availableWeight=0.0; previous=null;
        for (i=0; i<=lastEndpoint; i++) {
            for (j=0; j<=endpoints[i].lastOpen; j++) {
                currentInterval=endpoints[i].open[j];
                currentInterval.independentSetWeight=availableWeight+currentInterval.weight;
                currentInterval.independentSetPrevious=previous;
            }
            for (j=0; j<=endpoints[i].lastClosed; j++) {
                currentInterval=endpoints[i].closed[j];
                if (currentInterval.independentSetWeight>availableWeight) {
                    availableWeight=currentInterval.independentSetWeight;
                    previous=currentInterval;
                }
            }
        }
        maxWeight=availableWeight;
        if (maxWeight==0) {
            System.err.println("The optimal solution of sample "+sample+" "+(onHap1?"hap1":"hap2")+" has weight zero?! sampleWindowLast="+sampleWindowLast);        
            System.exit(1);
        }

        // Removing from the current haplotype every interval that is not in the
        // selected max-weight independent set.
        for (i=sampleWindowLast; i>=0; i--) {
            if (sampleWindow[i].independentSetWeight!=maxWeight) continue;
            currentInterval=sampleWindow[i];
            while (currentInterval!=null) {
                currentInterval.inIndependentSet=true;
                currentInterval=currentInterval.independentSetPrevious;
            }
            break;
        }
        error=true;
        if (onHap1) {
            markIndependentSetOverlaps(sample,true,false,sampleWindowLast);
            for (i=0; i<=sampleWindowLast; i++) {
                if (sampleWindow[i].inIndependentSet) { error=false; continue; }
                if (sampleWindow[i].overlapsIS_hap1) sampleWindow[i].genotypes[sample]=REMOVE_HAP1_D[sampleWindow[i].genotypes[sample]];
                else sampleWindow[i].genotypes[sample]=REMOVE_HAP1_0[sampleWindow[i].genotypes[sample]];
            }
        }
        else {
            markIndependentSetOverlaps(sample,false,true,sampleWindowLast);
            for (i=0; i<=sampleWindowLast; i++) {
                if (sampleWindow[i].inIndependentSet) { error=false; continue; }
                if (sampleWindow[i].overlapsIS_hap2) sampleWindow[i].genotypes[sample]=REMOVE_HAP2_D[sampleWindow[i].genotypes[sample]];
                sampleWindow[i].genotypes[sample]=REMOVE_HAP2_0[sampleWindow[i].genotypes[sample]];
            }
        }
        if (error) {
            System.err.println("No selected interval in the optimal solution of sample "+sample+" "+(onHap1?"hap1":"hap2")+"?!  maxWeight="+maxWeight+" sampleWindowLast="+sampleWindowLast);        
            System.exit(1);
        }
    }


    /**
     * Stores in $endpoints$ the sorted list of distinct first and last 
     * positions of all intervals in $sampleWindow$.
     *
     * Remark: the procedure assumes that $endpoints$ has already been 
     * initialized.
     */
    private static final void intervals2endpoints() {
        int i, j;
        Point tmpPoint;

        lastEndpoint=-1;
        for (i=0; i<=sampleWindowLast; i++) {
            lastEndpoint++;
            endpoints[lastEndpoint].clear();
            endpoints[lastEndpoint].pos=sampleWindow[i].first;
            endpoints[lastEndpoint].addOpen(sampleWindow[i]);
            lastEndpoint++;
            endpoints[lastEndpoint].clear();
            endpoints[lastEndpoint].pos=sampleWindow[i].last;
            endpoints[lastEndpoint].addClosed(sampleWindow[i]);
        }
        if (lastEndpoint>1) {
            // If there are <=2 points, they must be already sorted. 
            Arrays.sort(endpoints,0,lastEndpoint+1);
        }
        j=0;
        for (i=1; i<=lastEndpoint; i++) {
            if (endpoints[i].pos==endpoints[j].pos) endpoints[j].merge(endpoints[i]);
            else {
                j++;
                tmpPoint=endpoints[j];
                endpoints[j]=endpoints[i];
                endpoints[i]=tmpPoint;
            }
        }
        lastEndpoint=j;
    }
    
    
    /**
     * Sets the $overlapsIS_hap*$ flags of every element of $sampleWindow[0..
     * last]$. Such flags tell if an interval that is not in the selected max
     * weight independent set, overlaps with an element of the selected max
     * weight independent set on a given haplotype of $sample$.
     *
     * Remark: the procedure assumes that $window[0..last]$ is sorted by last
     * position.
     *
     * @param hap* TRUE=the procedure considers this haplotype.
     */
    private static final void markIndependentSetOverlaps(int sample, boolean hap1, boolean hap2, int last) {
        boolean onHap1, onHap2;
        int i, j;
        int first;
        
        if (hap1) {
            for (i=0; i<=last; i++) sampleWindow[i].overlapsIS_hap1=false;
        }
        if (hap2) {
            for (i=0; i<=last; i++) sampleWindow[i].overlapsIS_hap2=false;
        }
        for (i=last; i>=0; i--) {
            first=sampleWindow[i].first;
            onHap1=sampleWindow[i].onHap1(sample);
            onHap2=sampleWindow[i].onHap2(sample);
            for (j=i-1; j>=0; j--) {
                if (sampleWindow[j].last<first) break;
                if (sampleWindow[i].inIndependentSet && !sampleWindow[j].inIndependentSet) {
                    if (hap1 && onHap1 && sampleWindow[j].onHap1(sample)) sampleWindow[j].overlapsIS_hap1=true;
                    if (hap2 && onHap2 && sampleWindow[j].onHap2(sample)) sampleWindow[j].overlapsIS_hap2=true;
                }
                else if (sampleWindow[j].inIndependentSet && !sampleWindow[i].inIndependentSet) {
                    if (hap1 && onHap1 && sampleWindow[j].onHap1(sample)) sampleWindow[i].overlapsIS_hap1=true;
                    if (hap2 && onHap2 && sampleWindow[j].onHap2(sample)) sampleWindow[i].overlapsIS_hap2=true;
                }
            }
        }
    }
    
    
    
    
    // --------------------- WINDOW DRAWING PROCEDURES -------------------------
    
    /**
     * Drawing constants
     */
    private static final int PIXELS_PER_POS = 50;
    private static final int PIXELS_DELTA = 5;
    private static final int N_ROWS = 5*PIXELS_PER_POS;
    private static final Color COLOR_BACKGROUND = new Color(0x00FFFFFF);
    private static final Stroke STROKE_UNPHASED = new BasicStroke(1.0f,BasicStroke.CAP_ROUND,BasicStroke.JOIN_ROUND,10.0f,new float[] {10.0f},0.0f);
    private static final Stroke STROKE_PHASED = new BasicStroke(1.0f,BasicStroke.CAP_ROUND,BasicStroke.JOIN_ROUND,10.0f);
    private static final Stroke STROKE_BACKGROUND = new BasicStroke(0.0f);
    private static Color[] type2color;
    
    
    private static final void initType2color() {
        type2color = new Color[6];
        type2color[DEL] = new Color(0x006d9eeb);  // Blue
        type2color[INV] = new Color(0x0093c47d);  // Green
        type2color[DUP] = new Color(0x00ffd966);  // Yellow
        type2color[INS] = new Color(0x00cc0000);  // Red
        type2color[SNP] = new Color(0x00666666);  // Gray
        type2color[REPLACEMENT] = new Color(0x009900ff);  // Violet
    }
    
    
    private static final void drawWindow(int sample, int last, Graphics2D graphics, int nColumns, int first, Random random) {
        int i, x, y;
        int gt, width;
        
        graphics.setStroke(STROKE_BACKGROUND);
        graphics.setColor(COLOR_BACKGROUND);
        graphics.fillRect(0,0,nColumns,N_ROWS);
        for (i=0; i<=last; i++) {
            gt=window[i].genotypes[sample];
            if (gt==UNPHASED_00 || gt==PHASED_00 || gt==UNPHASED_0D || gt==UNPHASED_D0 || gt==UNPHASED_DD || gt==PHASED_0D || gt==PHASED_D0 || gt==PHASED_DD) continue;
            if (gt==UNPHASED_01 || gt==UNPHASED_10 || gt==UNPHASED_D1 || gt==UNPHASED_1D) graphics.setStroke(STROKE_UNPHASED);
            else graphics.setStroke(STROKE_PHASED);  // Also assigned to 1/1
            graphics.setColor(type2color[window[i].variantType]);
            x=(1+window[i].first-first)*PIXELS_PER_POS-PIXELS_DELTA+random.nextInt((PIXELS_DELTA)<<1);
            width=(window[i].last-window[i].first+1)*PIXELS_PER_POS;
            if (gt==UNPHASED_01 || gt==UNPHASED_10 || gt==UNPHASED_D1 || gt==UNPHASED_1D || gt==UNPHASED_11 || gt==PHASED_11) {
                y=PIXELS_PER_POS-PIXELS_DELTA+random.nextInt((PIXELS_DELTA)<<1);
                graphics.drawRect(x,y,width,3*PIXELS_PER_POS);
            }
            else if (gt==PHASED_10) {
                y=PIXELS_PER_POS-PIXELS_DELTA+random.nextInt((PIXELS_DELTA)<<1);
                graphics.drawRect(x,y,width,PIXELS_PER_POS);
            }
            else if (gt==PHASED_01) {
                y=PIXELS_PER_POS*3-PIXELS_DELTA+random.nextInt((PIXELS_DELTA)<<1);
                graphics.drawRect(x,y,width,PIXELS_PER_POS);
            }
        }
    }
    
        
    
    
    // -------------------------- DATA STRUCTURES  -----------------------------
    
    /**
     * A biallelic VCF record, represented as a weighted interval on the line.
     * The same object is intended to be reused by multiple VCF records.
     */
    private static class Interval implements Comparable {
        public static final int ORDER_INPUT = 0;
        public static final int ORDER_LAST_POS = 1;
        public static int order;
        
        public boolean isSV;  // True=SV, False=SNP/small indel.
        public int variantType;
        public int chr;
        public int first, last;  // One-based, inclusive.
        public int inputIndex;
        public double weight;  // Assumed to be non-negative
        public String[] vcfRecord;  // One cell per VCF column
        public int[] genotypes;  // Biallelic
        
        /**
         * Independent set variables
         */
        public boolean inIndependentSet;
        public double independentSetWeight;
        public Interval independentSetPrevious;
        public boolean overlapsIS_hap1, overlapsIS_hap2;
        
        
        public Interval() {
            isSV=false; variantType=-1; chr=-1; first=-1; last=-1; weight=0.0; 
            vcfRecord=null; genotypes=null; inputIndex=-1;
        }
        
        
        public final void clearIndependentSetVariables() {
            independentSetWeight=0.0; independentSetPrevious=null; inIndependentSet=false;
            overlapsIS_hap1=false; overlapsIS_hap2=false;
        }
        
        
        /**
         * @param record a line of a VCF file.
         */
        public final void init(String record) {
            int i;
            int pos, length, nSamples;
            String tmpString;
            
            this.vcfRecord=record.split(VCF_SEPARATOR);
            tmpString=getInfoField(this.vcfRecord[7],SVTYPE_STR);
            if (tmpString!=null && tmpString.length()!=0) {
                isSV=true;
                variantType=svType2Row(tmpString);
            }
            else {
                isSV=false;
                variantType=refAlt2Row(this.vcfRecord[3],this.vcfRecord[4]);
            }
            if (variantType==-1) {
                System.err.println("ERROR: this record encodes an unknown variant type: "+record);
                System.exit(1);
            }
            chr=string2contig(this.vcfRecord[0]);
            pos=Integer.parseInt(this.vcfRecord[1]);
            tmpString=getInfoField(this.vcfRecord[7],SVLEN_STR);
            if (tmpString!=null) {
                length=Integer.parseInt(tmpString);
                if (length<0) length=-length;
            }
            else if (variantType==REPLACEMENT) length=this.vcfRecord[3].length()-1;
            else length=Math.max(this.vcfRecord[3].length(),this.vcfRecord[4].length())-1;
            first=-1; last=-1;
            if (variantType==DEL || variantType==INV || variantType==DUP || variantType==REPLACEMENT) { first=pos+1; last=pos+length; }
            else if (variantType==INS) { first=pos; last=pos+1; }
            else if (variantType==SNP) { first=pos; last=pos; }
            nSamples=this.vcfRecord.length-9;
            if (genotypes==null || genotypes.length<nSamples) genotypes = new int[nSamples];
            for (i=0; i<nSamples; i++) genotypes[i]=gt2Id(this.vcfRecord[9+i]);
        }
        
        
        /**
         * Sets $weight$ to a value loaded from the INFO field or the SAMPLE
         * field, depending on global variable $WEIGHT_TAG_IN_SAMPLE_COLUMN$.
         *
         * Remark: if no value is found in the VCF record, $weight$ is 
         * arbitrarily set to one.
         *
         * Remark: weight is assumed to be non-negative throughout the program,
         * but everything should still work with negative weights.
         */
        private final void setWeight(int sample) {
            int i, j, p;
            int gtLength;
            String gt, value;
            
            value=null;
            if (WEIGHT_TAG_IN_SAMPLE_COLUMN) {
                p=vcfRecord[8].indexOf(WEIGHT_TAG);
                if (p>=0) {
                    j=0;
                    for (i=0; i<p; i++) {
                        if (vcfRecord[8].charAt(i)==GT_SEPARATOR) j++;
                    }
                    gt=vcfRecord[9+sample]; gtLength=gt.length();
                    for (i=0; i<gtLength; i++) {
                        if (gt.charAt(i)!=GT_SEPARATOR) continue;
                        j--;
                        if (j>0) continue;
                        p=gt.indexOf(GT_SEPARATOR,i+1);
                        value=p>=0?gt.substring(i+1,p):gt.substring(i+1);
                        break;
                    }
                }
            }
            else value=getInfoField(vcfRecord[7],WEIGHT_TAG);
            weight=value!=null?Double.parseDouble(value):1.0;
        }
        
        
        public final boolean isPresent(int sample) {
            return genotypes[sample]!=PHASED_00 && genotypes[sample]!=UNPHASED_00 && 
                   genotypes[sample]!=PHASED_DD && genotypes[sample]!=UNPHASED_DD && 
                   genotypes[sample]!=PHASED_D0 && genotypes[sample]!=UNPHASED_D0 && 
                   genotypes[sample]!=PHASED_0D && genotypes[sample]!=UNPHASED_0D;
        }
        
        
        public final boolean onHap1(int sample) {
            return genotypes[sample]==PHASED_10 || genotypes[sample]==PHASED_1D ||
                   genotypes[sample]==PHASED_11 || 
                   genotypes[sample]==UNPHASED_11;
        }
        
        
        public final boolean onHap2(int sample) {
            return genotypes[sample]==PHASED_01 || genotypes[sample]==PHASED_D1 ||
                   genotypes[sample]==PHASED_11 || 
                   genotypes[sample]==UNPHASED_11;
        }
        
        
        public final boolean isPhased(int sample) {
            return genotypes[sample]==PHASED_00 || genotypes[sample]==PHASED_01 || genotypes[sample]==PHASED_10 || genotypes[sample]==PHASED_11 ||
                   genotypes[sample]==PHASED_DD || genotypes[sample]==PHASED_D1 || genotypes[sample]==PHASED_1D || genotypes[sample]==PHASED_D0 || genotypes[sample]==PHASED_0D;
        }

        
        /**
         * @param nextInterval an interval that ends after $last$;
         * @param sample both this interval and $nextInterval$ are assumed to be
         * phased in $sample$;
         * @return TRUE iff this interval can precede $nextInterval$ in an 
         * independent set of $sample$.
         */
        public final boolean precedes(Interval nextInterval, int sample) {
            return N_GT_COLLISIONS[genotypes[sample]][nextInterval.genotypes[sample]]==0 || last<nextInterval.first;
        }
        
        
        /**
         * Prints the original VCF record, but using the GTs in $genotypes$.
         */
        public final void toVCF(BufferedWriter bw) throws IOException {
            int i, j, k, p, q;
            int gtIndex, gtLength;
            String gt;
            
            // Computing the index of the GT field in this record
            p=vcfRecord[8].indexOf(GT_STR);
            gtIndex=0;
            for (i=0; i<p; i++) {
                if (vcfRecord[8].charAt(i)==GT_SEPARATOR) gtIndex++;
            }
            
            // Printing the VCF record
            bw.write(vcfRecord[0]);
            for (i=1; i<=8; i++) { bw.write(VCF_SEPARATOR); bw.write(vcfRecord[i]); }
            for (i=0; i<genotypes.length; i++) {
                bw.write(VCF_SEPARATOR);
                gt=vcfRecord[9+i]; gtLength=gt.length();
                k=0; p=0;
                for (j=0; j<gtLength; j++) {
                    if (gt.charAt(j)==GT_SEPARATOR) {
                        k++;
                        if (k==gtIndex) { p=j+1; break; }
                    }
                }
                if (p>0) bw.write(gt.substring(0,p));
                bw.write(GT2STR[genotypes[i]]);
                q=gt.indexOf(GT_SEPARATOR,p);
                if (q>=0) bw.write(gt.substring(q));
            }
            bw.newLine();
        }
        
        
        public String toString() {
            return (isSV?"1":"0")+", "+variantType+", chr"+chr+"["+first+".."+last+"], weight="+weight;
        }
        
        
        public boolean equals(Object other) {
            final Interval otherInterval = (Interval)other;
            
            if (order==ORDER_LAST_POS) return last==otherInterval.last;
            else if (order==ORDER_INPUT) return inputIndex==otherInterval.inputIndex;
            return false;
        }
        
        
        public int compareTo(Object other) {
            final Interval otherInterval = (Interval)other;
            
            if (order==ORDER_LAST_POS) {
                if (last<otherInterval.last) return -1;
                else if (last>otherInterval.last) return 1;
                else return 0;
            }
            else if (order==ORDER_INPUT) {
                if (inputIndex<otherInterval.inputIndex) return -1;
                else if (inputIndex>otherInterval.inputIndex) return 1;
                else return 0;
            }
            return 0;
        }
        
    }
    
    
	private static final int svType2Row(String type) {
		if ( type.equalsIgnoreCase(DEL_STR) || 
			 type.equalsIgnoreCase(DEL_ME_STR)
		   ) return DEL;
		else if (type.equalsIgnoreCase(INV_STR)) return INV;
        else if ( type.equalsIgnoreCase(DUP_STR) ||
			      type.equalsIgnoreCase(DUP_TANDEM_STR) ||
				  type.equalsIgnoreCase(DUP_INT_STR) ||
                  type.equalsIgnoreCase(CNV_STR)
			    ) return DUP;
        else if ( type.equalsIgnoreCase(INS_STR) ||
                  type.equalsIgnoreCase(INS_ME_STR) ||
                  type.equalsIgnoreCase(INS_NOVEL_STR)
                ) return INS;
		else return -1;
	}


	private static final int refAlt2Row(String ref, String alt) {
        if (ref.length()==1) {
            if (alt.length()>1) return INS;
            else return SNP;
        }
        else {
            if (alt.length()==1) return DEL;
            else return REPLACEMENT;
        }
	}
    
    
    /**
     * An endpoint of an interval
     */
    private static class Point implements Comparable {
        private static final int CAPACITY = 2;  // Arbitrary
        public int pos;  // One-based
        public Interval[] open, closed;
        public int lastOpen, lastClosed;


        public Point() {
            open = new Interval[CAPACITY];
            closed = new Interval[CAPACITY];
            clear();
        }


        public final void clear() {
            pos=-1; lastOpen=-1; lastClosed=-1;
        }


        public final void addOpen(Interval interval) {
            lastOpen++;
            if (lastOpen==open.length) {
                Interval[] newArray = new Interval[open.length<<1];
                System.arraycopy(open,0,newArray,0,open.length);
                open=newArray;
            }
            open[lastOpen]=interval;
        }


        public final void addClosed(Interval interval) {
            lastClosed++;
            if (lastClosed==closed.length) {
                Interval[] newArray = new Interval[closed.length<<1];
                System.arraycopy(closed,0,newArray,0,closed.length);
                closed=newArray;
            }
            closed[lastClosed]=interval;
        }


        /**
         * Adds to $open,closed$ the intervals in the corresponding arrays of
         * $otherPoint$.
         */
        public final void merge(Point otherPoint) {
            int newLength;

            newLength=lastOpen+1+otherPoint.lastOpen+1;
            if (newLength>open.length) {
                Interval[] newArray = new Interval[newLength];
                System.arraycopy(open,0,newArray,0,lastOpen+1);
                open=newArray;
            }
            System.arraycopy(otherPoint.open,0,open,lastOpen+1,otherPoint.lastOpen+1);
            lastOpen=newLength-1;
            newLength=lastClosed+1+otherPoint.lastClosed+1;
            if (newLength>closed.length) {
                Interval[] newArray = new Interval[newLength];
                System.arraycopy(closed,0,newArray,0,lastClosed+1);
                closed=newArray;
            }
            System.arraycopy(otherPoint.closed,0,closed,lastClosed+1,otherPoint.lastClosed+1);
            lastClosed=newLength-1;
        }


        public boolean equals(Object other) {
            final Point otherPoint = (Point)other;
            return pos==otherPoint.pos;
        }


        public int compareTo(Object other) {
            final Point otherPoint = (Point)other;
            if (pos<otherPoint.pos) return -1;
            else if (pos>otherPoint.pos) return 1;
            return 0;
        }


        public String toString() {
            int i;
            String out;

            out="pos="+pos+"\n  open: ";
            for (i=0; i<=lastOpen; i++) out+=open[i].variantType+",";
            out+="\n  closed: ";
            for (i=0; i<=lastClosed; i++) out+=closed[i].variantType+",";
            return out+"\n";
        }
    }
    
    
    
    
    // ------------------------- BASIC VCF HANDLING ----------------------------
    
    private static final String VCF_SEPARATOR = "\t";
    private static final char GT_SEPARATOR = ':';
    private static final String GT_STR = "GT";
    public static final char COMMENT = '#';
    public static final String END_STR = "END";
    public static final String CIEND_STR = "CIEND";
    public static final String INFO_SEPARATOR = ";";
    public static final String SVTYPE_STR = "SVTYPE";
    public static final String SVLEN_STR = "SVLEN";
    public static final String CHR_STR = "chr";
    public static final int CHR_STR_LENGTH = CHR_STR.length();
    public static final String X_STR_PRIME = "X";
    public static final String Y_STR_PRIME = "Y";
    public static final String M_STR_PRIME = "M";
    public static final String MT_STR_PRIME = "MT";
    public static final String X_STR = CHR_STR+X_STR_PRIME;
    public static final String Y_STR = CHR_STR+Y_STR_PRIME;
    public static final String M_STR = CHR_STR+M_STR_PRIME;
    public static final String MT_STR = CHR_STR+MT_STR_PRIME;
    
    /**
     * SV types: labels used by callers.
     */
    public static final String DEL_STR = "DEL";
    public static final String DEL_ME_STR = "DEL:ME";
    public static final String DEL_INV_STR = "DEL/INV";
    public static final String INS_STR = "INS";
    public static final String INS_ME_STR = "INS:ME";
    public static final String INS_NOVEL_STR = "INS:NOVEL";
    public static final String DUP_STR = "DUP";
    public static final String DUP_TANDEM_STR = "DUP:TANDEM";
    public static final String DUP_INT_STR = "DUP:INT";
    public static final String INV_STR = "INV";
    public static final String INV_DUP_STR = "INVDUP";
    public static final String CNV_STR = "CNV";
    public static final String BND_STR = "BND";
    public static final String TRA_STR = "TRA";
    
    
    /**
     * @return NULL if $field$ does not occur in $str$.
     */
    private static final String getInfoField(String str, String field) {
        final int FIELD_LENGTH = field.length()+1;
        int p = str.indexOf(field+"=");
        if (p<0) return null;
        if (field.equalsIgnoreCase(END_STR)) {
            while (p>=2 && str.substring(p-2,p-2+CIEND_STR.length()).equalsIgnoreCase(CIEND_STR)) p=str.indexOf(field+"=",p+1);
            if (p<0) return null;
        }
        final int q = str.indexOf(INFO_SEPARATOR,p+FIELD_LENGTH);
        return str.substring(p+FIELD_LENGTH,q<0?str.length():q);
    }
    
    
    /**
     * @return one-based.
     */
    private static final int string2contig(String str) {
        if (str.length()>=CHR_STR_LENGTH && str.substring(0,CHR_STR_LENGTH).equalsIgnoreCase(CHR_STR)) {
            if (str.equalsIgnoreCase(X_STR)) return 23;
            else if (str.equalsIgnoreCase(Y_STR)) return 24;
            else if (str.equalsIgnoreCase(M_STR) || str.equalsIgnoreCase(MT_STR)) return 25;
            else if (str.length()<=CHR_STR_LENGTH+2) return Integer.parseInt(str.substring(CHR_STR_LENGTH));
            else return -1;
        }
        else {
            if (str.equalsIgnoreCase(X_STR_PRIME)) return 23;
            else if (str.equalsIgnoreCase(Y_STR_PRIME)) return 24;
            else if (str.equalsIgnoreCase(M_STR_PRIME) || str.equalsIgnoreCase(MT_STR_PRIME)) return 25;
            else if (str.length()<=2) return Integer.parseInt(str);
            else return -1;
        }
    }

}