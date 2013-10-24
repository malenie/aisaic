package edu.vt.cbil.bacom;

import java.io.*;
import java.util.*;
import static java.lang.Math.*;

import umontreal.iro.lecuyer.probdist.NormalDist;


/* Author: Bai Zhang   
 * Email:  baizhang@vt.edu  */


/**
 * <code>CopyNumberDetection</code> detects, segments, and stores the copy
 * number profile of a tumor sample. 
 *
 * @version      2010.0727       
 */
public class CopyNumberDetection {

    /**
     * estimated noise standard variance in the copy number data
     */ 
    public double globalStd;

    /**
     * chromosome name: 1, 2, ..., 22, MT, X, Y
     */
    public ArrayList<String> chromosome;

    /**
     * segment start position
     */
    public ArrayList<Integer> segStartPos;

    /**
     * segment end position
     */
    public ArrayList<Integer> segEndPos;

    /** 
     * Segment status. The values are the estimated segment means. 2 normal, 
     * <2 deletion, >2 amplification.
     */
    public ArrayList<Double> segStatus;

    /**
     * Chromosome start index. Indicates the start index of each chromosome.
     */
    public ArrayList<Integer> chromosomeStartIndex;

    /**
     * p-value threshold for copy number change detection
     */
    private final double PVALUE_THRESH = 0.05;

    /**
     * the minimum length of copy number segment to be considered
     */
    private final int MIN_SEG_LENGTH = 30;

    /**
     * the threshold to be considered normal, between <code>2 -
     * normalSegThreshold</code> and  <code>2 + normalSegThreshold</code>.
     */
    private double normalSegThreshold = 0.1;


    /**
     * A constructor
     * @param cnv     copy number data
     */
    public CopyNumberDetection(AffyCopyNumberData cnv) {

	globalStd = estimateStd(cnv.copyNumber);
	chromosome = new ArrayList<String>(); 
	segStartPos = new ArrayList<Integer>(); /* Begins from 0 */
	segEndPos = new ArrayList<Integer>();   /* Begins from 0 */
	segStatus = new ArrayList<Double>();
	chromosomeStartIndex = new ArrayList<Integer>();

	if (cnv.chipType.equals("GenomeWideSNP_6")) {
	    normalSegThreshold = 0.2;
	    
	} else {
	    normalSegThreshold = 0.1;
	}

	detection(cnv);
    }


    /**
     * A constructor.
     * @param cnv          copy number data
     * @param threshold    the threshold to be considered normal, between
     *                     <code>2 - normalSegThreshold</code> and <code>2 +
     *                     normalSegThreshold</code>. 
     */
    public CopyNumberDetection(AffyCopyNumberData cnv, double threshold) {

	globalStd = estimateStd(cnv.copyNumber);
	chromosome = new ArrayList<String>(); 
	segStartPos = new ArrayList<Integer>(); /* Begins from 0 */
	segEndPos = new ArrayList<Integer>();   /* Begins from 0 */
	segStatus = new ArrayList<Double>();
	chromosomeStartIndex = new ArrayList<Integer>();

	normalSegThreshold = threshold;

	detection(cnv);
    }


    /**
     * Detects DNA copy number changes.
     * @param cnv     copy number data
     */
    private void detection(AffyCopyNumberData cnv) {

	String chr = cnv.chrID[0];
	int chrStart = 0;
	int chrEnd = 0;
	NormalDist gaussian = new NormalDist();
	double[] cn;

	this.chromosomeStartIndex.add(0);
	    
	for (int i = 0; i < cnv.chrID.length; i ++) {


	    /* Perform segmentation chromosome by chromosome. */

	    if (!chr.equals(cnv.chrID[i])) {

		chrEnd = i - 1;

		cn = new double[chrEnd - chrStart + 1];

		for (int j = 0; j < cn.length; j ++) {
		    cn[j] = cnv.copyNumber[chrStart + j] - 2;    
		          // Place normal segment mean to zero
		}

		double threshold = abs(gaussian.inverseF(PVALUE_THRESH 
			       / ((double)cn.length * (double)(cn.length - 1)))); 
		// threshold of the z-score given the p-value cutoff
		// cn.length * (cn.length - 1) is used for correcting multiple
		// testing 

		segmentation(chr, cn, 0, 0, threshold);

		chrStart = i;
		chr = cnv.chrID[i];

		this.chromosomeStartIndex.add(chrStart);
	    }

	}


	cn = new double[cnv.copyNumber.length - chrStart];

	for (int j = 0; j < cn.length; j ++) {
	    cn[j] = cnv.copyNumber[chrStart + j] - 2;
	}

	double threshold = abs(gaussian.inverseF(PVALUE_THRESH 
                        / ((double)cn.length * (double)(cn.length - 1))));
	segmentation(chr, cn, 0, 0, threshold);


	/* Cleaning up the detection results */
	int segIndex = 0;
	for (int i = 0; i < chromosomeStartIndex.size(); i ++) {
	    int startLocation = chromosomeStartIndex.get(i);
	    while (cnv.chrID[startLocation].equals(chromosome.get(segIndex))) {
		cn = new double[segEndPos.get(segIndex) - segStartPos.get(segIndex) + 1];
		System.arraycopy(cnv.copyNumber, startLocation+segStartPos.get(segIndex), cn, 0, cn.length);
		double segMean = mean(cn);
		segStatus.add(abs(segMean - 2)<normalSegThreshold? 2: segMean);
		segIndex ++;
		if (segIndex == chromosome.size()) {
		    break;
		}
	    }

	}

    }


    
    /**
     * Performs copy number change segmentation
     * @param chr             chromosome
     * @param cn              copy number
     * @param startLocation   start location
     * @param initialStatus   initial status
     * @param threshold       threshold
     */
    private void segmentation(String chr, double[] cn, int startLocation, int initialStatus, double threshold) {
	
	double score = 0;
	int startPos = -1;
	int endPos = -1;
	double previousScore = 0;

	double meanCN = mean(cn);

	if (initialStatus != 0) {
	    for (int i = 0; i < cn.length; i ++) {
		cn[i] = cn[i] - meanCN;
	    }
	}


	int[] breakPoints = new int[2];
	fastSegmentation(cn, breakPoints, threshold);
	startPos = breakPoints[0];
	endPos = breakPoints[1];

	  
	/* Standard method to detect startPos and endPos*/
// 	double sigma = estimateStd(cn);
// 	double score_old = 0;
// 	int index = 0;
	
// 	for (int start = 0; start < (cn.length - 1); start ++) {

// 	    score = (cn[start] + cn[start+1]) / (sqrt(2) * sigma);
// 	    score_old = score;

// 	    for (int end = start + 2; end < cn.length; end ++) {
		
// 		score = (score_old *sqrt(end - start) * sigma + cn[end])
// 		    / (sqrt(end - start + 1) * sigma);
// 		score_old = score;

// 		if (abs(score) >= threshold && abs(score) >= abs(previousScore) && (end-start+1) >= MIN_SEG_LENGTH) {
// 		    startPos = start;
// 		    endPos = end;
// 		    previousScore = score;
// 		}

// 	    }

// 	}



	if (startPos >= 0) {

	    if (startPos >= MIN_SEG_LENGTH) {
		double[] seg1 = new double[startPos];
		System.arraycopy(cn, 0, seg1, 0, startPos);
		double meanSeg1 = mean(seg1);
		for (int i = 0; i < seg1.length; i ++) {
		    seg1[i] -= meanSeg1;
		}
		segmentation(chr, seg1, startLocation, initialStatus, threshold);
	    } else {

		if (startPos >=1) {
		    this.chromosome.add(chr);
		    this.segStartPos.add(startLocation);
		    this.segEndPos.add(startLocation + startPos - 1);
		}
	    }

	    if ((endPos - startPos + 1) >= MIN_SEG_LENGTH) {
		double[] seg2 = new double[endPos - startPos + 1];
		System.arraycopy(cn, startPos, seg2, 0, endPos - startPos + 1);
		double meanSeg2 = mean(seg2);
		for (int i = 0; i < seg2.length; i ++) {
		    seg2[i] -= meanSeg2;
		}
		segmentation(chr, seg2, startLocation+startPos, 
			     initialStatus == 0? (int)signum(meanSeg2) : initialStatus, threshold);
	    } else {
		this.chromosome.add(chr);
		this.segStartPos.add(startLocation + startPos);
		this.segEndPos.add(startLocation + endPos);
	    }

	    if ((cn.length - endPos - 1) >= MIN_SEG_LENGTH) {
		double[] seg3 = new double[cn.length - endPos - 1];
		System.arraycopy(cn, endPos+1, seg3, 0, cn.length - endPos -1);
		double meanSeg3 = mean(seg3);
		for (int i = 0; i < seg3.length; i ++) {
		    seg3[i] -= meanSeg3;
		}
		segmentation(chr, seg3, startLocation + endPos + 1, initialStatus, threshold);
	    } else {
		

		if ((endPos+1) <= (cn.length-1)) {
		    this.chromosome.add(chr);
		    this.segStartPos.add(startLocation + endPos + 1);
		    this.segEndPos.add(startLocation + cn.length - 1);
		}
	    }


	} else {
	    
		this.chromosome.add(chr);
		this.segStartPos.add(startLocation);;
		this.segEndPos.add(startLocation + cn.length - 1);
	}


    }



    /**
     * Performs fast copy number change segmentation
     * @param cn              copy number
     * @param breakPoints     break points
     * @param threshold       threshold
     */
    private void fastSegmentation(double cn[], int breakPoints[], double threshold) {

	double score = 0;
	int startPos = -1;
	int endPos = -1;
	double previousScore = 0;

	final int MERGED_SEG_LENGTH = 15;
	
	int len = (int) ceil((double) cn.length / (double) MERGED_SEG_LENGTH);
	
	double x[] = new double[len];

	for (int i = 0; i < (len-1); i ++) {
	    
	    x[i] = 0;
	    for (int j = i * MERGED_SEG_LENGTH; j < (i + 1) * MERGED_SEG_LENGTH; j ++) {
		x[i] += cn[j];
	    }
	    x[i] = x[i] / MERGED_SEG_LENGTH; 
	}
	x[len-1] = 0;
	for (int j = (len - 1) * MERGED_SEG_LENGTH; j < cn.length; j++) {
	    x[len - 1] += cn[j];
	}
	x[len-1] = x[len-1] / (cn.length - (len - 1) * MERGED_SEG_LENGTH);
	


	double sigma = estimateStd(x);
	double score_old = 0;
	int index = 0;

	NormalDist gaussian = new NormalDist();
	double threshold2 = abs(gaussian.inverseF(gaussian.cdf01(-threshold) 
						  * MERGED_SEG_LENGTH * MERGED_SEG_LENGTH));
	
	for (int start = 0; start < (x.length - 1); start ++) {

	    score = (x[start] + x[start+1]) / (sqrt(2) * sigma);
	    score_old = score;

	    for (int end = start + 2; end < x.length; end ++) {
		
		score = (score_old *sqrt(end - start) * sigma + x[end])
		    / (sqrt(end - start + 1) * sigma);
		score_old = score;

		if (abs(score) >= threshold2 && abs(score) >= abs(previousScore) 
		    && (end-start+1) >= (MIN_SEG_LENGTH / MERGED_SEG_LENGTH) ) {
		    startPos = start;
		    endPos = end;
		    previousScore = score;
		}

	    }

	}


	/* Now we need to test the exact positions of the start and the end of 
	   the segment in array cn. start in [a, b-1], end in [c, d-1], and the 
	   fixed part in [b, c-1].
	 */
	
	int a = max(0, (startPos - 1) * MERGED_SEG_LENGTH);
	int b = min(cn.length , (startPos + 1) * MERGED_SEG_LENGTH);
	int c = (endPos - 1) * MERGED_SEG_LENGTH;
	int d = min(cn.length, (endPos + 1) * MERGED_SEG_LENGTH);


	score = 0;
	startPos = -1;
	endPos = -1;
	previousScore = 0;

	double sumFixedPart = 0;

	for (int i = b; i < c; i++) {
	    sumFixedPart += cn[i];
	}

	sigma = estimateStd(cn);
	index = 0;
	double sum1 = sumFixedPart; 
	double sum2 = 0;
	
	len = (endPos  - startPos - 2) * MERGED_SEG_LENGTH;
	
	for (int start = b - 1; start >= a; start --) {
	    
	    sum1 = sum1 + cn[start];
	    
	    sum2 = sum1;

	    for (int end = max(start, c); end < d; end ++) {
		
		sum2 = sum2 + cn[end];
		score = sum2 / (sqrt(end - start + 1) * sigma);

		if (abs(score) >= threshold && abs(score) >= abs(previousScore) && (end-start+1) >= MIN_SEG_LENGTH) {
		    startPos = start;
		    endPos = end;
		    previousScore = score;
		}

	    }

	}


	breakPoints[0] = startPos;
	breakPoints[1] = endPos;
	

    }


    /**
     * Calculates the mean of <code>array x</code>.
     * @param x     an array
     * @return      the estimated mean of <code>x</code>
     */
    private double mean(double[] x) {

	double sum = 0;
	for (int i = 0; i < x.length; i ++) {
	    sum += x[i];
	}

	return (sum / x.length);
       
    }


    /**
     * Calculates the standard deviation of <code>array x</code>.
     * @param x     an array
     * @return      the estimated standard deviation of <code>x</code>
     */
    private double estimateStd(double[] cn) {

	double sum = 0;
	for (int i = 1; i < cn.length; i ++) {
	    sum += pow(cn[i] - cn[i-1], 2);
	}

	return sqrt(sum / (2*(cn.length - 2)));

    }


}






