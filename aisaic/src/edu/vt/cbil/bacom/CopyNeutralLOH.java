package edu.vt.cbil.bacom;

import java.io.*;
import java.util.*;
import static java.lang.Math.*;


import umontreal.iro.lecuyer.probdist.ChiSquareNoncentralDist;
import umontreal.iro.lecuyer.probdist.ChiSquareDist;
import umontreal.iro.lecuyer.probdist.NormalDist;


public class CopyNeutralLOH {

    public String[] chrID;
    public double[] alleleA;   // allale A signal of only normal loci with genotype AB
    public double[] alleleB;   // allele B signal of only normal loci with genotype AB
    public double[] copyNumber;  // copy number of only normal loci with genotype AB
    public int[] location;  // the physical location on the chromosome
    public boolean[] isGenotypeAB;
    public int[] indexOnChr;   // index on each chromosome
    public int[] normalSegIndex;  


    public ArrayList<String> chromosome;
    public ArrayList<Integer> cnnlohSegStartLocation;
    public ArrayList<Integer> cnnlohSegEndLocation;
    public ArrayList<Integer> cnnlohSegStartIndex;
    public ArrayList<Integer> cnnlohSegEndIndex;

    double globalStd;


    private final double PVALUE_THRESH = 0.001;

    private final int MIN_SEG_LENGTH = 30;


    public CopyNeutralLOH (AffyCopyNumberData cnv, CopyNumberDetection cnvDetection) {


	int numNormalLociAB = 0;
	int segIndex = 0;
	for (int i = 0; i < cnvDetection.chromosomeStartIndex.size(); i ++) {
	    int startLocation = cnvDetection.chromosomeStartIndex.get(i);
	    
	    while (cnv.chrID[startLocation].equals(cnvDetection.chromosome.get(segIndex))) {
		
		double tmpSegStatus = cnvDetection.segStatus.get(segIndex);

		if (tmpSegStatus == 2) {
		    for (int j = cnvDetection.segStartPos.get(segIndex); 
			 j <= cnvDetection.segEndPos.get(segIndex); j ++) {

			if (cnv.isGenotypeAB[startLocation + j]) {
			    numNormalLociAB ++;
			}

		    }
		}


		segIndex ++;
		if (segIndex == cnvDetection.chromosome.size()) {
		    break;
		}
	    }

	}



	chrID = new String[numNormalLociAB];
	alleleA = new double[numNormalLociAB];
	alleleB = new double[numNormalLociAB];
	copyNumber = new double[numNormalLociAB];
	location = new int[numNormalLociAB];
	isGenotypeAB = new boolean[numNormalLociAB];
	indexOnChr = new int[numNormalLociAB];
	normalSegIndex = new int[numNormalLociAB];

	segIndex = 0;
	int index = 0;
	for (int i = 0; i < cnvDetection.chromosomeStartIndex.size(); i ++) {
	    int startLocation = cnvDetection.chromosomeStartIndex.get(i);
	    
	    while (cnv.chrID[startLocation].equals(cnvDetection.chromosome.get(segIndex))) {
		
		double tmpSegStatus = cnvDetection.segStatus.get(segIndex);
		if (tmpSegStatus == 2) {
		    for (int j = cnvDetection.segStartPos.get(segIndex); 
			 j <= cnvDetection.segEndPos.get(segIndex); j ++) {

			if (cnv.isGenotypeAB[startLocation + j]) {
			    chrID[index] = cnv.chrID[startLocation + j];
			    alleleA[index] = cnv.alleleA[startLocation + j];
			    alleleB[index] = cnv.alleleB[startLocation + j];
			    copyNumber[index] = cnv.copyNumber[startLocation + j];
			    location[index] = cnv.location[startLocation + j];
			    isGenotypeAB[index] = cnv.isGenotypeAB[startLocation + j];
			    indexOnChr[index] = j;
			    normalSegIndex[index] = segIndex;

			    index ++;

			}

		    }
		}


		segIndex ++;
		if (segIndex == cnvDetection.chromosome.size()) {
		    break;
		}
	    }

	}


	globalStd = estimateStd(alleleA, alleleB);

	chromosome = new ArrayList<String>();
	cnnlohSegStartLocation = new ArrayList<Integer>();
	cnnlohSegEndLocation = new ArrayList<Integer>();
	cnnlohSegStartIndex = new ArrayList<Integer>();
	cnnlohSegEndIndex = new ArrayList<Integer>();

	copyNeutralLohDetection();

    }




    



    public void copyNeutralLohDetection() {

	int currSegIndex = normalSegIndex[0];
	int segStart = 0;
	int segEnd = 0;
	double threshold;
	double[] s;
	NormalDist gaussian = new NormalDist();
	
	for (int i = 0; i < chrID.length; i ++) {

	    if (normalSegIndex[i] != currSegIndex) {

		segEnd = i - 1;

		s = new double[segEnd - segStart + 1];

		for (int j = 0; j < s.length; j ++) {
		    s[j] = pow((alleleA[segStart + j] - alleleB[segStart + j]), 2);
		}
		
		if (s.length > MIN_SEG_LENGTH) {
		    threshold = abs(gaussian.inverseF(PVALUE_THRESH 
						      / (0.5 * (double)s.length * (double)(s.length - 1))));
		    cnnlohSegmentation(s, segStart, threshold);
		}

		segStart = i;
		currSegIndex = normalSegIndex[i];


	    }

	}


	s = new double[chrID.length - segStart];
	for (int j = 0; j < s.length; j ++) {
	    s[j] = pow((alleleA[segStart + j] - alleleB[segStart + j]), 2);
	}


	if (s.length > MIN_SEG_LENGTH) {
	    threshold = abs(gaussian.inverseF(PVALUE_THRESH / (0.5 * (double)s.length * (double)(s.length - 1))));
	    cnnlohSegmentation(s, segStart, threshold);
	}





    }


    private void cnnlohSegmentation(double[] s, int segStart, double threshold) {
	

	int[] breakPoints;

	breakPoints = regularSegmentation(s, threshold);
	//breakPoints = regularSegmentationChi2(s);
	
	int bp1 = breakPoints[0];
	int bp2 = breakPoints[1];

	if (bp1 != -1 && bp2 != -1) {

	    System.out.println("Chr = " + chrID[segStart] + "\n\n");

	    if (bp1 >= MIN_SEG_LENGTH) {
		double[] seg1 = new double[bp1];
		System.arraycopy(s, 0, seg1, 0, bp1);
		cnnlohSegmentation(seg1, segStart, threshold);
	    }

	    chromosome.add(chrID[segStart + bp1]);
	    cnnlohSegStartLocation.add(location[segStart + bp1]);
	    cnnlohSegEndLocation.add(location[segStart + bp2]);
	    cnnlohSegStartIndex.add(indexOnChr[segStart + bp1]);
	    cnnlohSegEndIndex.add(indexOnChr[segStart + bp2]);

	    if ((s.length - bp2 - 1) >= MIN_SEG_LENGTH) {
		double[] seg3 = new double[s.length - bp2 - 1];
		System.arraycopy(s, bp2+1, seg3, 0, s.length - bp2 - 1);
		cnnlohSegmentation(seg3, segStart, threshold);

	    }


	}




    }




    private int[] regularSegmentation(double[] s, double threshold) {

	/* Standard method to detect startPos and endPos*/

	int[] breakPoints = new int[2];
	double sigma = globalStd;
	double score = 0;
	double score_old = 0;
	int index = 0;
	double previousT = 0;

	
	int startPos = -1;
	int endPos = -1;


	for (int start = 0; start < (s.length - 1); start ++) {

	    score = (s[start] + s[start+1]) / (sigma*sigma);
	    score_old = score;

	    for (int end = start + 2; end < s.length; end ++) {
		
		score = (score_old  + s[end]/(sigma*sigma));
		score_old = score;

// 		double t = (score - (end-start+1)) / sqrt(2*(end-start+1));
		double t = (sqrt(score*2) - sqrt((end-start+1)*2 - 1));

		if (t >= threshold && t >= previousT && (end-start+1) >= MIN_SEG_LENGTH) {
		    startPos = start;
		    endPos = end;
		    previousT = t;
		}

	    }	
	}




	if (startPos != -1 && endPos != -1) {
	    NormalDist gaussian = new NormalDist();
	    System.out.println("startPos = " + startPos);
	    System.out.println("endPos = " + endPos);
	    System.out.println("previousT = " + previousT);
	    System.out.println("threshold = " + threshold);
	    System.out.println("pVal = " + (1-gaussian.cdf(previousT)));
	}


	breakPoints[0] = startPos;
	breakPoints[1] = endPos;

	return breakPoints;
	

    }



    /* The following code uses chi-square to detect the cnnloh segments.
       One problem is that it may not follow chi-square distribution exactly.
       The code above resorts to the asymptotic properties of chi-square / 
       chi distribution.*/     

    private int[] regularSegmentationChi2(double[] s) {

	/* Standard method to detect startPos and endPos*/

	int[] breakPoints = new int[2];
	double sigma = globalStd;
	double score = 0;
	double score_old = 0;
	int index = 0;
	double pVal = 1;
	double previousPVal = 1;
	double threshold = PVALUE_THRESH / ((s.length-1) * s.length / 2);
	
	int startPos = -1;
	int endPos = -1;


	for (int start = 0; start < (s.length - 1); start ++) {


	    score = (s[start] + s[start+1]) / (sigma*sigma);
	    score_old = score;

	    for (int end = start + 2; end < s.length; end ++) {
		
		score = (score_old  + s[end]/(sigma*sigma));
		score_old = score;
		
		ChiSquareDist chiSquare = new ChiSquareDist(end - start + 1);
		pVal = 1 - chiSquare.cdf(score);
		
		
		if (pVal <= threshold && pVal <= previousPVal && (end-start+1) >= MIN_SEG_LENGTH) {

		    startPos = start;
		    endPos = end;
		    previousPVal = pVal;
		}

	    }	

	}

	breakPoints[0] = startPos;
	breakPoints[1] = endPos;

	return breakPoints;
	

    }





    private double estimateStd(double[] x, double[] y) {

	double sum = 0;

	double mean = 0;
	for (int i = 0; i < x.length; i ++) {
	    mean += x[i]- y[i];
	}
	
	mean = mean / x.length;

	for (int i = 0; i < x.length; i ++) {
	    sum += pow(x[i] - y[i] - mean, 2);
	}

	return sqrt(sum / (x.length - 1));




    }
    



}





