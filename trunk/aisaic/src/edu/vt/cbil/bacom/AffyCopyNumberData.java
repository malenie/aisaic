package edu.vt.cbil.bacom;

import java.io.*;
import java.util.*;
import static java.lang.Math.*;

import affymetrix.fusion.cdf.FusionCDFData;
import affymetrix.fusion.cdf.FusionCDFProbeGroupInformation;
import affymetrix.fusion.cdf.FusionCDFProbeInformation;
import affymetrix.fusion.cdf.FusionCDFProbeSetInformation;
import affymetrix.fusion.cdf.FusionGeneChipProbeSetType;
import affymetrix.fusion.cel.FusionCELData;

import umontreal.iro.lecuyer.probdist.NormalDist;


/* Author: Bai Zhang   
 * Email:  baizhang@vt.edu  */


/**
 * <code>AffyCopyNumberData</code> extracts from CEL files, annotates, and
 * stores Affymetrix copy number data. An <code>AffyCopyNumberData</code> object
 * encapsulates the copy number data for further copy number analysis, such as
 * copy number detection, and normal tissue fraction estimation. 
 *
 * @version      2010.0803
 */
public class AffyCopyNumberData {

    /** 
     * an array of probe set IDs
     */
    public String[] probeSetID;

    /** 
     * an array of allele A signals
     */
    public double[] alleleA;

    /** 
     * an array allele B signals
     */
    public double[] alleleB;

    /** 
     * an array of copy number intensities from the normal sample
     */
    private double[] intensityNormal;

    /** 
     * an array of copy number intensities from the cancer sample
     */
    private double[] intensityTumor;

    /** 
     * an array of DNA copy numbers
     */
    public double[] copyNumber;

    /** 
     * an array of chromosome IDs
     */
    public String[] chrID;

    /** 
     * an array of physical locations on the corresponding chromosomes
     */
    public int[] location;

    /** 
     * an array of boolean values specifying whether it is genotype AB
     */
    public boolean[] isGenotypeAB;

    /** 
     * an array of boolean values specifying whether it is outlier
     */
    public boolean[] isOutlier;

    /** 
     * Affymetrix SNP chip type
     */
    public String chipType;

    /** 
     * an array of probe set types: SNP or CN
     */
    public int[] probeSetType;

    /** 
     * whether the files are successfully read
     */
    public boolean isSuccessfulOpenFile;

    /** 
     * whether the annotation files are successfully read
     */
    public boolean isSuccessfulOpenAnnotationFile;


    /** 
     * the ratio of average Y chromosome intensity and average X chromosome
     * intensity 
     */
    public double yxIntensityRatio;


    /** 
     * the number of annotated probe sets
     */
    private int numAnnotatedProbeSet;
    

    /**
     * A constructor for the class AffyCopyNumberData when Affymetrix Genomewide
     * 5.0, Affymetrix Genowide 6.0 are used.
     * 
     * @param normalCelFileName    Affymetrix CEL file of the paired normal
     *                             sample
     * @param tumorCelFileName     Affymetrix CEL file of the tumor sample
     * @param cdfFileName          Affymetrix library file (CDF file) for the
     *                             platform on which the samples are generated
     */
    public AffyCopyNumberData (String normalCelFileName, String tumorCelFileName, String cdfFileName, String annotationFileName) {

        /* Check if these files are from the same platform*/
	if (!checkChipType(normalCelFileName, tumorCelFileName, cdfFileName)) {
	    return;
	}

	try {
	    isSuccessfulOpenAnnotationFile = true;
	    readAnnotationFile(annotationFileName);
	} catch (IOException e){
	    isSuccessfulOpenAnnotationFile = false;
	}

	if (isSuccessfulOpenAnnotationFile && isSuccessfulOpenFile) {
	
	    readCel(normalCelFileName, tumorCelFileName, cdfFileName);

	    removeOutlier();

	    calculateYXratio();
	}

    }



    /**
     * A constructor for the class AffyCopyNumberData when Affymetrix 500K set
     * is used.
     * 
     * @param normalNspCelFileName    Affymetrix nsp CEL file of the paired
     *                                normal sample
     * @param normalStyCelFileName    Affymetrix sty CEL file of the paired
     *                                normal sample
     * @param tumorNspCelFileName     Affymetrix nsp CEL file of the paired
     *                                tumor sample
     * @param tumorStyCelFileName     Affymetrix sty CEL file of the paired
     *                                tumor sample
     * @param nspCdfFileName          Affymetrix library file (CDF file) for
     *                                Affymetrix 500K nsp chip
     * @param styCdfFileName          Affymetrix library file (CDF file) for
     *                                Affymetrix 500K sty chip
     */
    public AffyCopyNumberData(String normalNspCelFileName, 
			String normalStyCelFileName, String tumorNspCelFileName, 
			String tumorStyCelFileName, String nspCdfFileName, 
			String styCdfFileName, String annotationFileName) {

	if (!checkChipType500K(normalNspCelFileName, normalStyCelFileName, tumorNspCelFileName, 
			       tumorStyCelFileName, nspCdfFileName, styCdfFileName)) {
	    return;
	}

	try {
	    isSuccessfulOpenAnnotationFile = true;
	    readAnnotationFile(annotationFileName);
	} catch (IOException e){
	    isSuccessfulOpenAnnotationFile = false;
	}

	if (isSuccessfulOpenAnnotationFile && isSuccessfulOpenFile) {
	    readCel(normalNspCelFileName, tumorNspCelFileName, nspCdfFileName);
	    readCel(normalStyCelFileName, tumorStyCelFileName, styCdfFileName);
	
	    removeOutlier();

	    calculateYXratio();

	}
    }


    /**
     * Checks whether CEL files and CDF files are from the same platform.
     * 
     * @param normalCelFileName    Affymetrix CEL file of the paired normal
     *                             sample
     * @param tumorCelFileName     Affymetrix CEL file of the tumor sample
     * @param cdfFileName          Affymetrix library file (CDF file) for the
     *                             platform on which the samples are generated
     * @return                     <code>true</code> if the CEL files and CDF
     *                             file are from the same platform, otherwise
     *                             <code>false</code> 
     */
    private boolean checkChipType(String normalCelFileName, String tumorCelFileName, String cdfFileName) {

	// Check whether the normal CEL file, tumor CEL file, and the CDF file are consistent, and set
	// the variable chipType.

	FusionCELData celNormal;
	FusionCELData celTumor;
	FusionCDFData cdf;

	celNormal = new FusionCELData();
	celNormal.setFileName(normalCelFileName);
	if (celNormal.read() == false) {
	    isSuccessfulOpenFile = false;
	    return false;
	}

	celTumor = new FusionCELData();
	celTumor.setFileName(tumorCelFileName);
	if (celTumor.read() == false) {
	    isSuccessfulOpenFile = false;
	    return false;
	}
	
	cdf = new FusionCDFData();
	cdf.setFileName(cdfFileName);
	if (cdf.read() == false) {
	    isSuccessfulOpenFile = false;
	    return false;
	}


	if (!celNormal.getChipType().equals(cdf.getChipType()) ) {
	    isSuccessfulOpenFile = false;
	    return false;
	}

	if (!celTumor.getChipType().equals(cdf.getChipType()) ) {
	    isSuccessfulOpenFile = false;
	    return false;
	}

	chipType = celNormal.getChipType();
	isSuccessfulOpenFile = true;
	
	return true;
    }


    /**
     * Checks whether CEL files and CDF files are from the same platform
     * (Affymetrix 500K Array Set).
     * 
     * @param normalNspCelFileName    Affymetrix nsp CEL file of the paired
     *                                normal sample
     * @param normalStyCelFileName    Affymetrix sty CEL file of the paired
     *                                normal sample
     * @param tumorNspCelFileName     Affymetrix nsp CEL file of the paired
     *                                tumor sample
     * @param tumorStyCelFileName     Affymetrix sty CEL file of the paired
     *                                tumor sample
     * @param nspCdfFileName          Affymetrix library file (CDF file) for
     *                                Affymetrix 500K nsp chip
     * @param styCdfFileName          Affymetrix library file (CDF file) for
     *                                Affymetrix 500K sty chip
     * @return                        <code>true</code> if the CEL files and CDF
     *                                file are from the same platform, otherwise
     *                                <code>false</code> 
     */
    private boolean checkChipType500K(String normalNspCelFileName, String normalStyCelFileName, 
				      String tumorNspCelFileName, String tumorStyCelFileName, 
				      String nspCdfFileName, String styCdfFileName) {

	/* Check whether the normal CEL file, tumor CEL file, and the CDF file are from Affymetrix Human 
	 * Mapping 500K SNP chips (Nsp and Sty), and set the variable chipType to Mapping500K.
	 */

	if (!checkChipType(normalNspCelFileName, tumorNspCelFileName, nspCdfFileName) ) {
	    isSuccessfulOpenFile = false;
	    return false;
	}

	if (!chipType.equals("Mapping250K_Nsp") ) {
	    isSuccessfulOpenFile = false;
	    return false;
	}

	if ( !checkChipType(normalStyCelFileName, tumorStyCelFileName, styCdfFileName) ) {
	    isSuccessfulOpenFile = false;
	    return false;
	}

	if (!chipType.equals("Mapping250K_Sty") ) {
	    isSuccessfulOpenFile = false;
	    return false;
	}

	chipType = "Mapping500K";
	isSuccessfulOpenFile = true;

	return true;

    }


    /**
     * Reads the annotation file. Annotation files are located under the folder
     * annotation of the software directory.
     * @throws IOException
     */
    private void readAnnotationFile(String annotationFileName) throws IOException {

	// Read and process annotation files.
//	
//	String annotationFileName = "annotation/" + chipType + "_Annotation.csv";
//        

	numAnnotatedProbeSet = 0;

	BufferedReader annotationFile = null;
	
        try {
            annotationFile = new BufferedReader(new FileReader(annotationFileName));

	    String rowAnnotation;

	    while ((rowAnnotation = annotationFile.readLine()) != null) {
		numAnnotatedProbeSet ++;
            }

	} catch(Exception e) {
	    isSuccessfulOpenAnnotationFile = false;
	}finally {
	    if (annotationFile != null) {
                annotationFile.close();
            }
        }

	probeSetID = new String[numAnnotatedProbeSet];
	alleleA = new double[numAnnotatedProbeSet];
	alleleB = new double[numAnnotatedProbeSet];

	intensityNormal = new double[numAnnotatedProbeSet];
	intensityTumor = new double[numAnnotatedProbeSet];

	copyNumber = new double[numAnnotatedProbeSet];
	chrID = new String[numAnnotatedProbeSet];
	location = new int[numAnnotatedProbeSet];
	isGenotypeAB = new boolean[numAnnotatedProbeSet];
	probeSetType = new int[numAnnotatedProbeSet];
	isOutlier = new boolean[numAnnotatedProbeSet];

	for (int i = 0; i < isOutlier.length; i ++) {
	    isOutlier[i] = true;
	}

	annotationFile = null;
	
        try {
            annotationFile = new BufferedReader(new FileReader(annotationFileName));

	    String rowAnnotation;
	    int rowNum = 0;

	    while ((rowAnnotation = annotationFile.readLine()) != null) {
		String[] rowSplit = null;
		rowSplit = rowAnnotation.split(",");
		
		probeSetID[rowNum] = rowSplit[0];
		chrID[rowNum] = rowSplit[1];
		location[rowNum] = Integer.parseInt(rowSplit[2]);
		rowNum++;
		
            }
	} catch(Exception e) {
	    isSuccessfulOpenAnnotationFile = false;
	}
        finally {
	    if (annotationFile != null) {
                annotationFile.close();
            }
        }


    }


    /**
     * Reads CEL files. <code>readCel</code> reads the intensity signals from
     * the tumor sample and the normal sample. The copy number at a particular
     * locus is calculated using the ratio of the tumor intensity at that locus
     * and the correspoding normal intensity times 2.
     * 
     * @param normalCelFileName    Affymetrix CEL file of the paired normal
     *                             sample
     * @param tumorCelFileName     Affymetrix CEL file of the tumor sample
     * @param cdfFileName          Affymetrix library file (CDF file) for the
     *                             platform on which the samples are generated
     */
    private void readCel(String normalCelFileName, String tumorCelFileName, String cdfFileName) {
	
	FusionCELData celNormal;
	FusionCELData celTumor;
	FusionCDFData cdf;
	
	celNormal = new FusionCELData();
	celNormal.setFileName(normalCelFileName);
	if (celNormal.read()== false) {
	    System.out.println("Failed to read the CEL file.");
	    return;
	}

	celTumor = new FusionCELData();
	celTumor.setFileName(tumorCelFileName);
	if (celTumor.read()== false) {
	    System.out.println("Failed to read the CEL file.");
	    return;
	}
	
	cdf = new FusionCDFData();
	cdf.setFileName(cdfFileName);
	if (cdf.read() == false) {
	    System.out.println("Failed to read the CDF file.");
	    return;
	}

	int nsets = cdf.getHeader().getNumProbeSets();

	ProbeSetIntensityData[] probeSetDataNormal = new ProbeSetIntensityData[nsets];
	ProbeSetIntensityData[] probeSetDataTumor = new ProbeSetIntensityData[nsets];
// SNP array structure:
//      1.  Each probeset contains several groups (the number of groups varies among different probeset)
//      2.  Each group contains several cells (the number of cells also varies among different groups)
 	for (int iset = 0; iset < nsets; iset++) {

	    String probeSetName = cdf.getProbeSetName(iset); //get the probeset name
	    FusionCDFProbeSetInformation set = new FusionCDFProbeSetInformation();
	    cdf.getProbeSetInformation(iset, set);
	    int ngroups = set.getNumGroups();

	    int numPmANormal = 0;  //Pm : perfect match
	    int numPmBNormal = 0;

	    int numPmATumor = 0;
	    int numPmBTumor = 0;


	    int numMmANormal = 0; // Mm: Mis-Match
	    int numMmBNormal = 0;

	    int numMmATumor = 0;
	    int numMmBTumor = 0;


	    probeSetDataNormal[iset] = new ProbeSetIntensityData();
	    probeSetDataTumor[iset] = new ProbeSetIntensityData();


	    probeSetDataNormal[iset].probeSetType =set.getProbeSetType();
	    probeSetDataTumor[iset].probeSetType =set.getProbeSetType();

	    probeSetDataNormal[iset].probeSetID = cdf.getProbeSetName(iset)+"";
	    probeSetDataTumor[iset].probeSetID = cdf.getProbeSetName(iset)+"";

	    for (int igroup = 0; igroup < ngroups; igroup++) {

		FusionCDFProbeGroupInformation group = new FusionCDFProbeGroupInformation();
		set.getGroup(igroup, group);
		int ncells = group.getNumCells();

		for (int icell = 0; icell < ncells; icell++) {

		    FusionCDFProbeInformation probe = new FusionCDFProbeInformation();
		    group.getCell(icell, probe);

		    try {

			char pBase = probe.getPBase();
			char tBase = probe.getTBase();
//                       only If the match is perfect, the intensity of this cell contributes
			if ((((pBase + tBase) == 213) || ((pBase + tBase) == 202))) {
//            Perfect match is the match that with pBase:tBase = a:t or c:g
			    if ((igroup % 2) == 0) {

				if (!celNormal.isOutlier(probe.getX(), probe.getY())) {
				    probeSetDataNormal[iset].pmA += celNormal.getIntensity(probe.getX(), probe.getY());
				    numPmANormal ++;
				}

				if (!celTumor.isOutlier(probe.getX(), probe.getY())) {
				    probeSetDataTumor[iset].pmA += celTumor.getIntensity(probe.getX(), probe.getY());
				    numPmATumor ++;
				}

				
			    } else {

				if (!celNormal.isOutlier(probe.getX(), probe.getY())) {
				    probeSetDataNormal[iset].pmB += celNormal.getIntensity(probe.getX(), probe.getY());
				    numPmBNormal ++;
				}

				if (!celTumor.isOutlier(probe.getX(), probe.getY())) {
				    probeSetDataTumor[iset].pmB += celTumor.getIntensity(probe.getX(), probe.getY());
				    numPmBTumor ++;
				}
			    }


			} 
		    } catch(Exception e) {
		    }
		}

	    }
//using the average of intensity of the perfect match cells as the intensity for certain probeset

	    if (numPmANormal != 0) {
		probeSetDataNormal[iset].pmA = probeSetDataNormal[iset].pmA / numPmANormal;
	    }
	    if (numPmBNormal != 0) {
		probeSetDataNormal[iset].pmB = probeSetDataNormal[iset].pmB / numPmBNormal;
	    }
	    if (numPmATumor != 0) {
		probeSetDataTumor[iset].pmA = probeSetDataTumor[iset].pmA / numPmATumor;
	    }
	    if (numPmBTumor != 0) {
		probeSetDataTumor[iset].pmB = probeSetDataTumor[iset].pmB / numPmBTumor;
	    }


	}

	try{
	    genotypeCalling(probeSetDataNormal);  //  to decide whether certain probe is AB type
	} catch (Exception e) {}

	Arrays.sort(probeSetDataNormal);
	Arrays.sort(probeSetDataTumor);

	ArrayList<Integer> indexFoundList = new ArrayList<Integer>(probeSetDataNormal.length);
	ArrayList<Integer> indexFoundListSNP = new ArrayList<Integer>(probeSetDataNormal.length);

	for (int i = 0; i < numAnnotatedProbeSet; i++) {
	    int indexFound = Arrays.binarySearch(probeSetDataNormal, new ProbeSetIntensityData(probeSetID[i]));

	    if (indexFound >=0) {

		isGenotypeAB[i] = probeSetDataNormal[indexFound].isGenotypeAB;
		probeSetType[i] = probeSetDataNormal[indexFound].probeSetType;

		copyNumber[i] = 2 * (probeSetDataTumor[indexFound].pmA + probeSetDataTumor[indexFound].pmB) 
		    / (probeSetDataNormal[indexFound].pmA + probeSetDataNormal[indexFound].pmB) ;

		intensityNormal[i] = probeSetDataNormal[indexFound].pmA + probeSetDataNormal[indexFound].pmB;

		intensityTumor[i] = probeSetDataTumor[indexFound].pmA + probeSetDataTumor[indexFound].pmB;

//                For 250k affymetrix chip, all the probes are SNP probes
		if (probeSetType[i] == FusionGeneChipProbeSetType.GenotypingProbeSetType) {
		    alleleA[i] = probeSetDataTumor[indexFound].pmA / probeSetDataNormal[indexFound].pmA;
		    alleleB[i] = probeSetDataTumor[indexFound].pmB / probeSetDataNormal[indexFound].pmB;

		    if(alleleA[i] > 1E-10  && alleleA[i] <= 30 && alleleB[i] >= 1E-10 && alleleB[i] <= 30
		       && copyNumber[i] > 1E-10 && copyNumber[i] <= 30) {
			isOutlier[i] = false;
			indexFoundList.add(i);
			indexFoundListSNP.add(i);
		    }
		} 
//              For SNP6.0, there are half SNP probes and half CN probes
 		if (probeSetType[i] == FusionGeneChipProbeSetType.CopyNumberProbeSetType) {
		    if(copyNumber[i] > 1E-10 && copyNumber[i] <= 30) {
			isOutlier[i] = false;
			indexFoundList.add(i);
		    }
		}  

	    }
	}

	/* Median filtering on the intensity data before calculating the copy numbers*/
	// Filters f1 = new Filters(intensityNormal);
	// f1.medianFilter(3);
	// Filters f2 = new Filters(intensityTumor);
	// f2.medianFilter(3);



	// for (int j = 0; j < indexFoundList.size(); j ++) {

	//     int i = indexFoundList.get(j);

	//     copyNumber[i] = intensityTumor[i] / intensityNormal[i] * 2;

	//     if (!(copyNumber[i] > 1E-10 && copyNumber[i] <=30)) {
	// 	copyNumber[i] = 2.0;
	// 	isOutlier[i] = true;
	//     }

	// }

	globalNormalization(indexFoundList, indexFoundListSNP);

    }


    /**
     * Makes genotype AB calls. <code>genotypeCalling</code> determines which
     * loci are genotype AB (heterozygous) based on A and B allelesl signal of
     * the normal sample.
     * 
     * @param p    probe set intensity data
     */
    private void genotypeCalling(ProbeSetIntensityData[] p) throws IOException {

	int numGenotypingProbeSet = 0;
	double meanPmA = 0;
	double meanPmB = 0;

	for (int i = 0; i < p.length; i ++) {
	    if (p[i].probeSetType == FusionGeneChipProbeSetType.GenotypingProbeSetType) {
		meanPmA += (p[i].pmA == 0 ? 0:log2(p[i].pmA));
		meanPmB += (p[i].pmB == 0 ? 0:log2(p[i].pmB));
		numGenotypingProbeSet ++;

	    }
	}

	meanPmA = meanPmA / numGenotypingProbeSet;
	meanPmB = meanPmB / numGenotypingProbeSet;

	double k = meanPmB / meanPmA;


	class DistanceToLine implements Comparable<DistanceToLine> {
	    public Double dist;
	    public int index;
	    public boolean equals(DistanceToLine d) {
		return (dist == d.dist);
	    }
	    public int compareTo(DistanceToLine d) {
		return dist.compareTo(d.dist);
	    }
	}


	DistanceToLine[] s = new DistanceToLine[numGenotypingProbeSet];
	int j = 0;
	for (int i = 0; i < p.length; i ++) {
	    if (p[i].probeSetType == FusionGeneChipProbeSetType.GenotypingProbeSetType) {
		s[j] = new DistanceToLine();
		s[j].dist = Math.abs(log2(p[i].pmB) - k * log2(p[i].pmA));
		s[j].index = i;
		j ++;
	    }
	}

	Arrays.sort(s);
//  Add for testing
//         BufferedWriter genotypeCallingFile = null;
//try {
//   genotypeCallingFile = new BufferedWriter(new FileWriter("genoTypeCalling_s_dist.csv"), 30*1024*1024);
//
//   for (int i = 0; i < numGenotypingProbeSet; i ++) {
//       String entry = s[i].dist + "," + s[i].index;
//       genotypeCallingFile.write(entry.toCharArray());
//        genotypeCallingFile.newLine();
//   }
//
//   } finally {
//       if (genotypeCallingFile != null) {
//           genotypeCallingFile.close();
//       }
//}
////
	for (int i = 0; i < Math.ceil(numGenotypingProbeSet * 0.20); i ++)
        {
	    
	    p[s[i].index].isGenotypeAB = true;

	}


// 	BufferedWriter genotypeCallingFile = null;
// 	try {
//             genotypeCallingFile = new BufferedWriter(new FileWriter("genoTypeCallingData.csv"), 30*1024*1024);

// 	    for (int i = 0; i < p.length; i ++) {
// 		if (p[i].probeSetType == FusionGeneChipProbeSetType.GenotypingProbeSetType) {
// 		    String entry = (p[i].pmA == 0 ? 0:log2(p[i].pmA)) + "," + (p[i].pmB == 0 ? 0:log2(p[i].pmB))
// 			+ "," + (p[i].isGenotypeAB? 1:0);
// 		    genotypeCallingFile.write(entry.toCharArray());
// 		    genotypeCallingFile.newLine();
// 		}
// 	    }

// 	} finally {
// 	    if (genotypeCallingFile != null) {
//                 genotypeCallingFile.close();
//             }
//         }




    }


    /** 
     * Removes outliers. This method is deprecated. It removes loci whose
     * intensities are outliers. The resultant copy number vectors are of
     * different lengths.
     */ 
    private void removeOutlier_old() {
	
	int numValidProbeSets = 0;
	for (int i = 0; i < probeSetID.length; i ++) {
	    if (!isOutlier[i]) {
		numValidProbeSets ++;
	    }
	}
	String[] probeSetID_new = new String[numValidProbeSets];
	double[] alleleA_new = new double[numValidProbeSets];
	double[] alleleB_new = new double[numValidProbeSets];

	double[] intensityNormal_new = new double[numValidProbeSets];
	double[] intensityTumor_new = new double[numValidProbeSets];

	double[] copyNumber_new = new double[numValidProbeSets];
	String[] chrID_new = new String[numValidProbeSets];
	int[] location_new = new int[numValidProbeSets];
	boolean[] isGenotypeAB_new = new boolean[numValidProbeSets];
	boolean[] isOutlier_new = new boolean[numValidProbeSets];
	int[] probeSetType_new = new int[numValidProbeSets];

	int j = 0;
	for (int i = 0; i < probeSetID.length; i ++) {
	    if (!isOutlier[i]) {
		probeSetID_new[j] = probeSetID[i];
		alleleA_new[j] = alleleA[i];
		alleleB_new[j] = alleleB[i];

		intensityNormal_new[j] = intensityNormal[i];
		intensityTumor_new[j] = intensityTumor[i];

		copyNumber_new[j] = copyNumber[i];
		chrID_new[j] = chrID[i];
		location_new[j] = location[i];
		isGenotypeAB_new[j] = isGenotypeAB[i];
		isOutlier_new[j] = isOutlier[i];
		probeSetType_new[j] = probeSetType[i];
		
		j++;
		
	    }
	}

	probeSetID = probeSetID_new;
	alleleA = alleleA_new;
	alleleB = alleleB_new;

	intensityNormal = intensityNormal_new;
	intensityTumor = intensityTumor_new;

	copyNumber = copyNumber_new;
	chrID = chrID_new;
	location = location_new;
	isGenotypeAB = isGenotypeAB_new;
	isOutlier = isOutlier_new;
	probeSetType = probeSetType_new;
	
	/* Median filtering of the copy number signals */
	// Filters f = new Filters(copyNumber);
	// f.medianFilter(3);
		
    }


    /** 
     * Imputes outliers by averaging adjacent probe sets.
     */ 
    private void removeOutlier() {
	
	for (int i = 0; i < probeSetID.length; i ++) {

	    if (isOutlier[i]) {
		isGenotypeAB[i] = false;
		double cn_new = 0;

		double num = 0;
		for (int j = max(-3, -i); j < 0; j ++) {
                    if (!isOutlier[i + j]) {
                        cn_new = cn_new + copyNumber[i + j];
                        num++;
                    }
		}
		for (int j = 1; j <= min(3, probeSetID.length - 1 - i); j ++) {
                    if (!isOutlier[i + j]) {
                        cn_new = cn_new + copyNumber[i + j];
                        num++;
                    }
		}
	        
                if (num >= 1 ) {
                    copyNumber[i] = cn_new / num;
                } else {
                    copyNumber[i] = 2;
                }
                alleleA[i] = 0;
                alleleB[i] = 0;
	    }
	    
	}

    }


    /**
     * Calculates the ratio of average Chromosome Y intensity and average
     * Chromosome X intensity in order to determine the gender of the subject.
     */ 
    private void calculateYXratio() {
	double chrXIntensity = 0;
	double chrYIntensity = 0;
	double numChrXProbeSet = 0;
	double numChrYProbeSet = 0;
	for (int i = 0; i < probeSetID.length; i ++) {
	    if (chrID[i].equals("X") && !Double.isInfinite(intensityNormal[i]) && !Double.isNaN(intensityNormal[i])) {
		chrXIntensity = chrXIntensity + intensityNormal[i];
		numChrXProbeSet ++;
	    }

	    if (chrID[i].equals("Y") && !Double.isInfinite(intensityNormal[i]) && !Double.isNaN(intensityNormal[i])) {
		chrYIntensity = chrYIntensity + intensityNormal[i];
		numChrYProbeSet ++;
	    }

	}

	yxIntensityRatio = (chrYIntensity / numChrYProbeSet) / (chrXIntensity / numChrXProbeSet);

    }


    /**
     * Normalizes the copy number signal, allele A signal and allele B signal
     * globally. 
     */
    public void globalNormalization(ArrayList<Integer> index, ArrayList<Integer> indexSNP) {


	double globalMean = 0;
	double globalMeanA = 0;
	double globalMeanB = 0;
	double[] copyNumberOneChip = new double[index.size()];
 	double[] alleleAOneChip = new double[indexSNP.size()];
 	double[] alleleBOneChip = new double[indexSNP.size()];

	for (int i = 0; i < index.size(); i++) {
	    copyNumberOneChip[i] = copyNumber[index.get(i)];
	}

	for (int i = 0; i < indexSNP.size(); i++) {
 	    alleleAOneChip[i] = alleleA[indexSNP.get(i)];
 	    alleleBOneChip[i] = alleleB[indexSNP.get(i)];
	}

	
	globalMean = robustMu(copyNumberOneChip);
 	globalMeanA = robustMu(alleleAOneChip);
 	globalMeanB = robustMu(alleleBOneChip);

	for (int i = 0; i < index.size(); i ++) {
	    copyNumber[index.get(i)] = copyNumber[index.get(i)] / globalMean * 2;
 	    alleleA[index.get(i)] = alleleA[index.get(i)] / globalMeanA;
 	    alleleB[index.get(i)] = alleleB[index.get(i)] / globalMeanB;
	}
    


    }



    /**
     * Normalizes the copy number signal, allele A signal and allele B signal
     * chromosome by chromosome. 
     */
    public void localNormalization() {


	double chrMean = 0;
	double[] chrCopyNumber = null;


	String chr = chrID[0] + "";
	int chrStartPos = 0;
	for (int i = 0; i < copyNumber.length; i ++) {
	    
	    if(!chr.equals(chrID[i])) {

		chrCopyNumber = new double[i - chrStartPos];
		System.arraycopy(copyNumber, chrStartPos, chrCopyNumber, 0, i - chrStartPos);
		chrMean = robustMu(chrCopyNumber);
// Local normalization is done only when the difference is small, otherwise, it may result in some artifial
		if (abs(2 - chrMean) < 0.05) {
		    for (int j = chrStartPos; j < i; j ++) {
		    		    
			copyNumber[j] = copyNumber[j] / chrMean * 2;

		    }
		}

		chrStartPos = i;
		chr = chrID[i] + "";

	    }

	}

	chrCopyNumber = new double[copyNumber.length - chrStartPos];
	System.arraycopy(copyNumber, chrStartPos, chrCopyNumber, 0, copyNumber.length - chrStartPos);
	chrMean = robustMu(chrCopyNumber);

	if (abs(2 - chrMean) < 0.075) {
	    for (int j = chrStartPos; j < copyNumber.length; j ++) {
		    		    
		copyNumber[j] = copyNumber[j] / chrMean * 2;

	    }
	}



// 	double chrMean = 0;
// 	double[] chrCopyNumber = null;


// 	String chr = chrID[0] + "";
// 	int chrStartPos = 0;
// 	int chrLength = 0;
// 	for (int i = 0; i < copyNumber.length; i ++) {
	    
// 	    if(!chr.equals(chrID[i])) {

// 		chrCopyNumber = new double[chrLength];
// 		int k = 0;
// 		for (int j = chrStartPos; j < i; j ++) {
// 		    if (!isOutlier[j]) {
// 			chrCopyNumber[k] = copyNumber[j];
// 			k ++;
// 		    }
// 		}

// 		chrMean = robustMu(chrCopyNumber);

// 		if (abs(2 - chrMean) < 0.1) {
// 		    for (int j = chrStartPos; j < i; j ++) {
		    		    
// 			copyNumber[j] = copyNumber[j] / chrMean * 2;

// 		    }
// 		}

// 		chrStartPos = i;
// 		chrLength = 0;
// 		chr = chrID[i] + "";

// 	    }

// 	    if (!isOutlier[i]) {
// 		chrLength ++;
// 	    }
	    

// 	}

// 	chrCopyNumber = new double[chrLength];
// 	int k = 0;
// 	for (int j = chrStartPos; j < copyNumber.length; j ++) {
// 	    if (!isOutlier[j]) {
// 		chrCopyNumber[k] = copyNumber[j];
// 		k ++;
// 	    }
// 	}

// 	chrMean = robustMu(chrCopyNumber);

// 	if (abs(2 - chrMean) < 0.1) {
// 	    for (int j = chrStartPos; j < copyNumber.length; j ++) {
		    		    
// 		copyNumber[j] = copyNumber[j] / chrMean * 2;

// 	    }
// 	}

    
    }


    /** 
     * Estimates the mean of array <code>y</code>. <code>robustMu</code> uses a
     * Gaussian model to remove outliers iteratively.
     * @param y    an array
     * @return     the estimated mean of <code>y</code>
     */ 
    private double robustMu(double[] y) {


	ArrayList<Double> x = new ArrayList<Double>(y.length);

	for (int i = 0; i < y.length; i ++) {
	    x.add(y[i]);
	}


	double sigma = 0;
	NormalDist gaussian = new NormalDist();
	double cutoff = 0;
	double cutoff1 = 0;
	double cutoff2 = 0;

	int xSize = x.size();
	
	sigma = std(x);
	double areaRemoved = 0.5/xSize;

	cutoff = abs(gaussian.inverseF(areaRemoved));
	
	cutoff1 =  -cutoff * sigma + mean(x);
	cutoff2 = cutoff * sigma + mean(x);

	while (true) {


	    Iterator<Double> iter = x.iterator();
	    
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

	    
	    sigma = std(x);
	    sigma = sigma / sqrt( (1 - 2*areaRemoved - 2*cutoff*exp(-pow(cutoff,2)/2)/sqrt(2*PI)) / (1 - 2*areaRemoved) );

	    xSize = x.size();
	    areaRemoved = 0.5/xSize;

	    cutoff = abs(gaussian.inverseF(areaRemoved));

	    cutoff1 =  -cutoff * sigma + mean(x);
	    cutoff2 = cutoff * sigma + mean(x);

	    
 	}

	return mean(x);
    }


    /**
     * Calculates the mean of <code>ArrayList</code> x.
     * @param x     an <code>ArrayList</code>
     * @return      the estimated mean of this <code>ArrayList</code>
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
     * Calculates the standard deviation of <code>ArrayList</code> x.
     * @param x     an <code>ArrayList</code>
     * @return      the estimated standard deviation of this
     *              <code>ArrayList</code> 
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
     * Calculates the standard deviation of <code>ArrayList</code> x, given its mean.
     * @param x     an <code>ArrayList</code>
     * @return      the estimated standard deviation of this
     *              <code>ArrayList</code> 
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
     * Calculates the binary logarithm.
     * @param x    a <code>double</code> value
     * @return     log2(x)
     */
    private double log2(double x) {
	return Math.log(x)/ Math.log(2);
    }

}




/** 
 * <code>ProbeSetIntensityData</code> stores intensity data (pmA, pmB) from cel
 * files. <code>ProbeSetIntensityData</code> can be sorted according to
 * probeSetID. 
 */
class ProbeSetIntensityData implements Comparable<ProbeSetIntensityData> {
    public String probeSetID; 
    public double pmA;      /* perfect match intensity A-allele */
    public double pmB;      /* perfect match intensity B-allele */
    public double mmA;
    public double mmB;
    public int probeSetType;    /* probeSetType 2: SNP  5: copy number */
    public boolean isGenotypeAB;    /* it is true if the genotype of this probe set is AB */

    public ProbeSetIntensityData () {
	this.probeSetID = null;
	this.pmA = 0;
	this.pmB = 0;
	this.mmA = 0;
	this.mmB = 0;
	this.probeSetType = 0;
	this.isGenotypeAB = false;
    }

    public ProbeSetIntensityData (String s) {
	this.probeSetID = s;
	this.pmA = 0;
	this.pmB = 0;
	this.probeSetType = 0;
	this.isGenotypeAB = false;
    }

    public boolean equals(ProbeSetIntensityData d) {
	return probeSetID.equals(d.probeSetID);
    }

    public boolean equals(String s) {
	return probeSetID.equals(s);
    }

    public int compareTo(ProbeSetIntensityData d) {
	return probeSetID.compareTo(d.probeSetID);
    }

	    
}

