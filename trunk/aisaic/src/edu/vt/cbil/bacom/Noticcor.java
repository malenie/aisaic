package edu.vt.cbil.bacom;

import java.io.*;
import java.util.*;
import static java.lang.Math.*;

import umontreal.iro.lecuyer.probdist.ChiSquareNoncentralDist;
import umontreal.iro.lecuyer.probdist.ChiSquareDist;
import umontreal.iro.lecuyer.probdist.NormalDist;



/**
 * <code>Noticcor</code> differentiates DNA copy number deletion types, and
 * corrects normal tissue contamination in the tumor sample. Noticcor is short
 * for NOrmal TIssue Contamination CORrection.
 *
 * @version      2010.0727       
 */
public class Noticcor {


    /** 
     * chromosome information
     */
    String[] chrInfo;

    /** 
     * A allele copy number
     */
    double[] cnXa;

    /** 
     * B allele copy number
     */
    double[] cnXb;

    /** 
     * DNA total copy number
     */
    double[] copyNumber;

    /** 
     * DNA total copy number
     */
    double[] cncLabel;

    /** 
     * whether it is genotype AB
     */
    boolean[] isGenotypeAB;

    /**
     * the number of all loci
     */
    int numLoci;

    /**
     * the estimated normal tissue fraction
     */
    public double alpha;

    /**
     * deletion segment chromosome ID
     */
    public ArrayList<String> deletionSegChrID;

    /**
     * deletion segment start location
     */
    public ArrayList<Integer> deletionSegStartLocation;

    /**
     * deletion segment end location
     */
    public ArrayList<Integer> deletionSegEndLocation;

    /**
     * the posterior probabilities that deletion segments are hemi-deletion
     */
    public ArrayList<Double> deletionSegProbHemiDeletion;

    /**
     * the posterior probabilities that deletion segments are homo-deletion
     */
    public ArrayList<Double> deletionSegProbHomoDeletion;

    /**
     * the deletion segment based normal tissue fraction
     */
    public ArrayList<Double> deletionSegNormalTissueFraction;
    

    /** 
     * A constructor.
     * @param cnv             an <code>AffyCopyNumberData</code> object that
     *                        stores the Affymetrix copy number data.
     * @param cnvDetection    an <code>CopyNumberDetection</code> object that
     *                        stores the copy number segmentation results.
     */
    public Noticcor(AffyCopyNumberData cnv, CopyNumberDetection cnvDetection) {

	this.numLoci = cnv.chrID.length;
	this.chrInfo = cnv.chrID;
	this.cnXa = cnv.alleleA;
	this.cnXb = cnv.alleleB;
	this.copyNumber = cnv.copyNumber;
	this.isGenotypeAB = cnv.isGenotypeAB;

	cncLabel = new double[numLoci];	
	for (int i = 0; i < numLoci; i ++) {
	    cncLabel[i] = 0;
	}

	int segIndex = 0;
	for (int i = 0; i < cnvDetection.chromosomeStartIndex.size(); i ++) {
	    int startLocation = cnvDetection.chromosomeStartIndex.get(i);
	    
	    while (cnv.chrID[startLocation].equals(cnvDetection.chromosome.get(segIndex))) {
		
		double tmpSegStatus = cnvDetection.segStatus.get(segIndex);

		for (int j = cnvDetection.segStartPos.get(segIndex); 
		     j <= cnvDetection.segEndPos.get(segIndex); j ++) {
		    cncLabel[startLocation + j] = tmpSegStatus - 2;
		}

		segIndex ++;
		if (segIndex == cnvDetection.chromosome.size()) {
		    break;
		}
	    }

	}


	deletionSegChrID = new ArrayList<String>();
	deletionSegStartLocation = new ArrayList<Integer>();
	deletionSegEndLocation = new ArrayList<Integer>();
	deletionSegProbHemiDeletion = new ArrayList<Double>();
	deletionSegProbHomoDeletion = new ArrayList<Double>();
	deletionSegNormalTissueFraction = new ArrayList<Double>();


	normalization();
	
    }


    /** 
     * Performs global normalization and chromosome-based normalization for the
     * copy number data.
     */
    private void normalization() {

	double avgNormalCN = 0;
	double avgNormalA = 0;
	double avgNormalB = 0;
	double numNormalLoci = 0;
	double numNormalLociAB = 0;

	/* Global normaliztion */

	for (int i = 0; i < numLoci; i ++) {

	    if (abs(cncLabel[i]) <= Double.MIN_VALUE) {
		avgNormalCN += copyNumber[i];
		numNormalLoci ++;
		if (isGenotypeAB[i]) {
		    avgNormalA += cnXa[i];
		    avgNormalB += cnXb[i];
		    numNormalLociAB ++;
		}
	    }

	}

	avgNormalCN = avgNormalCN / numNormalLoci;
	avgNormalA = avgNormalA / numNormalLociAB;
	avgNormalB = avgNormalB / numNormalLociAB;


// 	System.out.println("Global avgNormalCN = " + avgNormalCN);
// 	System.out.println("Global avgNormalA = " + avgNormalA);
// 	System.out.println("Global avgNormalB = " + avgNormalB);
// 	System.out.println("Global numNormalLoci = " + numNormalLoci);
// 	System.out.println("Global numNormalLociAB = " + numNormalLociAB);

	for (int i = 0; i < numLoci; i ++) {

	    copyNumber[i] = copyNumber[i] / avgNormalCN * 2;
	    cnXa[i] = cnXa[i] / avgNormalA;
	    cnXb[i] = cnXb[i] / avgNormalB;

	}


	/* Local normaliztion */
	
	int index = 0;
	while (index < numLoci) {
	    String chr = chrInfo[index];
	    int chrStartPos = index;
	    avgNormalCN = 0;
	    avgNormalA = 0;
	    avgNormalB = 0;
	    numNormalLoci = 0;
	    numNormalLociAB = 0;
	    if (abs(cncLabel[index]) <= Double.MIN_VALUE) {
		avgNormalCN += copyNumber[index];
		numNormalLoci ++;
		if (isGenotypeAB[index]) {
		    avgNormalA += cnXa[index];
		    avgNormalB += cnXb[index];
		    numNormalLociAB ++;
		}

	    }

	    index ++;		

	    while (chr.equals(chrInfo[index]) && index < numLoci) {

		if (abs(cncLabel[index]) <= Double.MIN_VALUE) {
		    avgNormalCN += copyNumber[index];
		    numNormalLoci ++;
		    if (isGenotypeAB[index]) {
			avgNormalA += cnXa[index];
			avgNormalB += cnXb[index];
			numNormalLociAB ++;
		    }
		}
		index ++;		
		if (index >= numLoci) {
		    break;
		}
	    }

	    if (numNormalLociAB > 500) {   /* if the number of normal loci with genotyp AB is greater
					      than 500, then perform normalization.*/

		avgNormalCN = avgNormalCN / numNormalLoci;
		avgNormalA = avgNormalA / numNormalLociAB;
		avgNormalB = avgNormalB / numNormalLociAB;

// 		System.out.println("Chr" + chrInfo[index-1] + " avgNormalCN = " + avgNormalCN);
// 		System.out.println("Chr" + chrInfo[index-1] + " avgNormalA = " + avgNormalA);
// 		System.out.println("Chr" + chrInfo[index-1] + " avgNormalB = " + avgNormalB);
// 		System.out.println("Chromosome numNormalLoci = " + numNormalLoci);
// 		System.out.println("Chromosome numNormalLociAB = " + numNormalLociAB);

		for (int i = chrStartPos; i < index; i++) {
		    copyNumber[i] = copyNumber[i] / avgNormalCN * 2;
		    cnXa[i] = cnXa[i] / avgNormalA;
		    cnXb[i] = cnXb[i] / avgNormalB;
		}
	    }

	    
	}

    }



    /**
     * Corrects normal tissue contamination in DNA copy number data.
     */ 
    public void correctNTC() {
	int startPos = 0;
	int endPos = -1;
	int m = 0;
	ArrayList<Double> segCN = new ArrayList<Double>(10000);
	ArrayList<Double> segXa = new ArrayList<Double>(10000);
	ArrayList<Double> segXb = new ArrayList<Double>(10000);

	double muAB;
	Double[] posteriorProb;
	double rhoAB = 0;
	String atChr = null;
	
	rhoAB = robustRho();

	alpha = 0;
	double sumSegSize = 0;

	for (int i = 1; i < numLoci; i++) {

	    segCN.clear();
	    segXa.clear();
	    segXb.clear();

	    if (cncLabel[i-1] < 0 && 
		(abs(cncLabel[i] - cncLabel[i-1]) > 0.000001 || i == (numLoci-1))) {
		

		if (i == (numLoci-1)) {
		    endPos = i;
		} else {
		    endPos = i - 1;
		}

		for (int j = startPos; j <= endPos; j++) {
		    segCN.add(copyNumber[j]);
		    if (isGenotypeAB[j]) {
			segXa.add(cnXa[j]);
			segXb.add(cnXb[j]);
		    }
		}
		

		if (segXa.size() >= 5 && (!chrInfo[endPos].equals("MT"))
		    && (!chrInfo[endPos].equals("X")) && (!chrInfo[endPos].equals("Y")) ) {
		    muAB = robustMu(segCN);
		    posteriorProb = bayesianPosterior(segCN, segXa, segXb, muAB, rhoAB, 0.01);
		    atChr = chrInfo[endPos];

		    if (!Double.isNaN(posteriorProb[2])) {
			alpha = alpha + (1 - posteriorProb[2]) * posteriorProb[3];
			sumSegSize = sumSegSize + posteriorProb[3];
		    }

		    
		    deletionSegChrID.add(atChr);
		    deletionSegStartLocation.add(startPos);
		    deletionSegEndLocation.add(endPos);
		    deletionSegProbHemiDeletion.add(posteriorProb[0]);
		    deletionSegProbHomoDeletion.add(posteriorProb[1]);
		    deletionSegNormalTissueFraction.add(1 - posteriorProb[2]);

		}
		
	    } 
	    
	    if (cncLabel[i] < 0 && 
		(abs(cncLabel[i] - cncLabel[i-1]) > Double.MIN_VALUE)) {

		startPos = i;
	    }
	}

        alpha = alpha / sumSegSize;


    }


    /**
     * Estimates the copy number segment mean. The robust estimation is achieved
     * by removing outliers iteratively.
     */
    public double robustMu(ArrayList<Double> segCN) {

	ArrayList<Double> sumCNaCNb = new ArrayList<Double>(segCN.size());
	double sigma = 0;
	NormalDist gaussian = new NormalDist();
	double segSize = 0;
	double cutoff = 0;
	double cutoff1 = 0;
	double cutoff2 = 0;

	for (int i = 0; i< segCN.size(); i++) {  // Remove outliers before compute mu
	    if (segCN.get(i) >= -50 && segCN.get(i) <= 50) {
		sumCNaCNb.add(segCN.get(i));
	    }
	}

	segSize = sumCNaCNb.size();
	
	sigma = std(sumCNaCNb);
	double areaRemoved = 0.5/segSize;

	cutoff = abs(gaussian.inverseF(areaRemoved));
	
	cutoff1 =  -cutoff * sigma + mean(sumCNaCNb);
	cutoff2 = cutoff * sigma + mean(sumCNaCNb);

	while (true && segSize >=10) {


	    Iterator<Double> iter = sumCNaCNb.iterator();
	    
	    int numRemoved = 0;
	    while (iter.hasNext()) {

		double cn = iter.next();
		if (cn < cutoff1 || cn > cutoff2) {
		    numRemoved++;
		    iter.remove();
		}
	    }

	    if (numRemoved == 0) {
		break;
	    }

	    
	    sigma = std(sumCNaCNb);
	    sigma = sigma / sqrt( (1 - 2*areaRemoved - 2*cutoff*exp(-pow(cutoff,2)/2)/sqrt(2*PI)) / (1 - 2*areaRemoved) );

	    segSize = sumCNaCNb.size();
	    areaRemoved = 0.5/segSize;

	    cutoff = abs(gaussian.inverseF(areaRemoved));

	    cutoff1 =  -cutoff * sigma + mean(sumCNaCNb);
	    cutoff2 = cutoff * sigma + mean(sumCNaCNb);

	    
 	}

	return mean(sumCNaCNb);
    }



    /**
     * Estimates the correlation between A-allele and B-allele signals for
     * genotype AB loci.The robust estimation is achieved by removing outliers
     * iteratively. 
     */
    private double robustRho() {

	ArrayList<Double> sumCNaCNb = new ArrayList<Double>(cnXa.length);
	ArrayList<Double> sXa = new ArrayList<Double>(cnXa.length);
	ArrayList<Double> sXb = new ArrayList<Double>(cnXa.length);
	double sigma = 0;
	NormalDist gaussian = new NormalDist();
	double segSize = 0;
	double cutoff = 0;
	double cutoff1 = 0;
	double cutoff2 = 0;

	for (int i = 0; i< cnXa.length; i++) {  // Remove outliers before compute mu
	    if (copyNumber[i] >= -50 && copyNumber[i] <= 50 && isGenotypeAB[i] && abs(cncLabel[i]) <= 0.001) {
		sumCNaCNb.add(copyNumber[i]);
		sXa.add(cnXa[i]);
		sXb.add(cnXb[i]);
	    }
	}

	segSize = (double) sumCNaCNb.size();
	
	sigma = std(sumCNaCNb);
	double areaRemoved = 0.5/segSize;

	cutoff = abs(gaussian.inverseF(areaRemoved));
	
	cutoff1 =  -cutoff * sigma + mean(sumCNaCNb);
	cutoff2 = cutoff * sigma + mean(sumCNaCNb);

	while (true && segSize >= 10) {


	    Iterator<Double> iter = sumCNaCNb.iterator();
	    Iterator<Double> iterXa = sXa.iterator();
	    Iterator<Double> iterXb = sXb.iterator();

	    int numRemoved = 0;
	    while (iter.hasNext()) {

		double cn = iter.next();
		iterXa.next();
		iterXb.next();
		if (cn < cutoff1 || cn > cutoff2) {
		    numRemoved++;
		    iter.remove();
		    iterXa.remove();
		    iterXb.remove();
		}
	    }

	    if (numRemoved == 0) {
		break;
	    }

	    
	    sigma = std(sumCNaCNb);
	    sigma = sigma / sqrt( (1 - 2*areaRemoved - 2*cutoff*exp(-pow(cutoff,2)/2)/sqrt(2*PI)) / (1 - 2*areaRemoved) );

	    segSize = sumCNaCNb.size();
	    areaRemoved = 0.5/segSize;

	    cutoff = abs(gaussian.inverseF(areaRemoved));

	    cutoff1 =  -cutoff * sigma + mean(sumCNaCNb);
	    cutoff2 = cutoff * sigma + mean(sumCNaCNb);
	    
 	}

	double corr = 0;
	double meanXa = mean(sXa);
	double meanXb = mean(sXb);

	for (int i = 0; i < sXa.size(); i++) {
	    corr = corr + (sXa.get(i) - meanXa) * (sXb.get(i) - meanXb);
	}

	corr = (corr / (sXa.size() - 1)) / (std(sXa) * std(sXb));

	return corr;
    }





    /** 
     * Calculates the mean of <code>ArrayList x<code>.
     * @param x    an <code>ArrayList<code>
     * @return     the estimated mean of <code>x</code>
     */
    private double mean(ArrayList<Double> x) {
	double sum = 0;
	double length = x.size();
	for (int i = 0; i < x.size(); i++) {
	    sum = sum + x.get(i);
	}
	return (sum / length);
    }


    /** 
     * Calculates the standard deviation of <code>ArrayList x<code>.
     * @param x    an <code>ArrayList<code>
     * @return     the estimated standard deviation of <code>x</code>
     */
    private double std(ArrayList<Double> x) {
	double sum = 0;
	double meanX = mean(x);
	double length = x.size();
	for (int i = 0; i < x.size(); i++) {
	    sum = sum + pow((x.get(i)-meanX), 2);
	}
	return (sqrt(sum / (length-1)));
    }


    /**
     * Calculates the standard deviation of <code>ArrayList x</code>, given its
     * mean. 
     * @param x       an <code>ArrayList<code>
     * @param meanX   average of <code>ArrayList x</code>
     * @return        the estimated standard deviation of <code>x</code>
     */
    private double std(ArrayList<Double> x, double meanX) {
	double sum = 0;
	double length = x.size();
	for (int i = 0; i < x.size(); i++) {
	    sum = sum + pow((x.get(i)-meanX), 2);
	}
	return (sqrt(sum / (length)));
    }


    /**
     * Calculates the Bayesian posterior probability of the status of the
     * deletion segment.
     * @param segCN       an <code>ArrayList<code> of copy number signals in the
     *                    segment
     * @param segXa       an <code>ArrayList</code> of A allele signals in the
     *                    segment
     * @param segXb       an <code>ArrayList</code> of B allele signals in the
     *                    segment 
     * @param muAB        segment mean
     * @param rhoAB       correlation between allele A and allele B
     * @param pTheshold   p-value threshold
     * @return            an <code>array</code> of pHemi (probability of
     *                    hemi-deletion), pHomo (probability of homo-deletion),
     *                    fractionTumor (estimated fraction of tumor), segSize
     *                    (segment size)
     */
    public Double[] bayesianPosterior(ArrayList<Double> segCN, ArrayList<Double> segXa, ArrayList<Double> segXb, 
				      double muAB, double rhoAB, double pThreshold) {

	ArrayList<Double> sumXaXb = new ArrayList<Double>(segXa.size());
	ArrayList<Double> sXa = new ArrayList<Double>(segXa.size());
	ArrayList<Double> sXb = new ArrayList<Double>(segXb.size());
	int segSize = segXa.size();

        for (int i = 0; i < segSize; i++) {
	    sumXaXb.add(segCN.get(i));
	    sXa.add(segXa.get(i));
	    sXb.add(segXb.get(i));
	}
	
	// *******************************
	// To remove outliers in the data
	// *******************************
	double sigma = std(sumXaXb, muAB);
	double areaRemoved = 0.5/segSize;
	NormalDist gaussian = new NormalDist();
	double cutoff = 0;
	double cutoff1 = 0;
	double cutoff2 = 0;

	cutoff = abs(gaussian.inverseF(areaRemoved));
	
	cutoff1 =  -cutoff * sigma + muAB;
	cutoff2 = cutoff * sigma + muAB;

	while (true && segSize >= 10) {

	    Iterator<Double> iter = sumXaXb.iterator();
	    Iterator<Double> iterXa = sXa.iterator();
	    Iterator<Double> iterXb = sXb.iterator();

	    int numRemoved = 0;
	    while (iter.hasNext()) {

		double cn = iter.next();
		iterXa.next();
		iterXb.next();
		if (cn < cutoff1 || cn > cutoff2) {
		
		    numRemoved++;
		    iter.remove();
		    iterXa.remove();
		    iterXb.remove();

		}
	    }

	    if (numRemoved == 0) {
		break;
	    }

	    sigma = std(sumXaXb, muAB);
	    sigma = sigma / sqrt( (1 - 2*areaRemoved - 2*cutoff*exp(-pow(cutoff,2)/2)/sqrt(2*PI)) / (1 - 2*areaRemoved) );


	    segSize = sumXaXb.size();

	    areaRemoved = 0.5/segSize;

	    cutoff = abs(gaussian.inverseF(areaRemoved));
	    cutoff1 =  -cutoff * sigma + muAB;
	    cutoff2 = cutoff * sigma + muAB;


 	}

	// *******************************

	segSize = sumXaXb.size();
	areaRemoved = 0.5/segSize;

	sigma = std(sumXaXb, muAB);
 	
	sigma = sigma / sqrt( (1 - 2*areaRemoved - 2*cutoff*exp(-pow(cutoff,2)/2)/sqrt(2*PI)) / (1 - 2*areaRemoved) );
	sigma = sigma * sqrt((1 - 1*rhoAB) / (1 + 1*rhoAB)); // corrected sigma
	
	double lambda = (2 - muAB) / sigma;
	lambda = pow(lambda, 2);


	// ***********************************
	// New method, using the whole segment
	// ***********************************

	double score = 0;
	for (int i = 0; i < segSize; i++) {
	    score = score + pow(sXa.get(i) - sXb.get(i), 2) / pow(sigma, 2);
	}

	ChiSquareDist chiSquare = new ChiSquareDist(segSize-1);

	ChiSquareNoncentralDist nx2 = new ChiSquareNoncentralDist(segSize - 1, lambda * segSize);
	Double L1 = nx2.density(score);
	Double L2 = chiSquare.density(score);
	
// 	System.out.println("mu = " + muAB);
// 	System.out.println("sigma = " + sigma);
//  	System.out.println("Score = " + score);
//  	System.out.println("segSize = " + segSize);
//  	System.out.println("lambda * segSize = " + (lambda * segSize));
// 	System.out.println("L1 = " + L1);
// 	System.out.println("L2 = " + L2);

	Double pHemi = new Double(0);
	Double pHomo = new Double(0);

	if (Double.isNaN(L1) && !Double.isNaN(L2)) {
	    pHemi = 0.0;
	    pHomo = 1.0;

	} 
	
	if (!Double.isNaN(L1) && Double.isNaN(L2)) {
	    pHemi = 1.0;
	    pHomo = 0.0;

	} 
	
	if (!Double.isNaN(L1) && !Double.isNaN(L2)) {
	    if (L1 == 0 && L2 != 0) {
		pHemi = 0.0;
		pHomo = 1.0;
	    } 
	    if (L2 == 0 && L1 != 0) {
		pHemi = 1.0;
		pHomo = 0.0;
	    }
	    if (L2 != 0 && L1 != 0) {
		pHemi = L1 / (L1 + L2);
		pHomo = L2 / (L1 + L2);
	    }
	    if (L2 == 0 && L1 == 0) {
		pHemi = 1.0;
		pHomo = 0.0;
	    }

	    
	}


	if (Double.isNaN(L1) && Double.isNaN(L2)) {
	    pHemi = 1.0;
	    pHomo = 0.0;

	}


	Double fractionTumor = Double.NaN;
	if (pHomo <= pThreshold) {   // Hemi-deletion
	    fractionTumor = 2 - muAB;   
	}
	if (pHemi <= pThreshold) {  // Homo-deletion
	    fractionTumor = 1 - muAB / 2;
	}

	Double result[] = new Double[4];
	result[0] = pHemi;
	result[1] = pHomo;
	result[2] = fractionTumor;
	result[3] = new Double(segSize);
	return result;
    }




}



