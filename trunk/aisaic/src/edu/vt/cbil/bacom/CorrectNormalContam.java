package edu.vt.cbil.bacom;

import edu.vt.cbil.util.BacomParameters;
import edu.vt.cbil.util.GeneDeletionInfo;
import edu.vt.cbil.util.Gene;
import edu.vt.cbil.util.Logger;
import java.io.*;
import java.util.concurrent.*;
import java.util.ArrayList;


public class CorrectNormalContam extends RecursiveAction {
   
    private int startSample;
    
    private int endSample;
        
    private BacomParameters bacomPara;
    
    
    public CorrectNormalContam(BacomParameters bacomPara, int start, int end)
    {
        this.bacomPara = bacomPara;
        this.startSample = start;
        this.endSample = end;
    }
    
    protected void compute()
    {
        /*
         * if parallel computing is wanted, these codes should be commented out  
         */
        for (int i = startSample; i<=endSample; i++)
        {
            BACOM(i, bacomPara.getNormalSamples().get(i), bacomPara.getTumorSamples().get(i), bacomPara.getInputCdfFile(),
                   bacomPara.getAnnotationFile(), bacomPara.getResultFolder(),  bacomPara.getSampleIds().get(i),  bacomPara.getGeneInterested());
        }
        /*
         * if parallel computing is wanted, these codes should be commented out  
         */
        //////////////////////////////////////////////////////////////////////////////////////
        /*
         * if serial computing is wanted, these codes should be commented out  
         */
//        if (startSample > endSample)
//        {
//            System.err.println("Error!");
//        }
//        if (startSample==endSample)
//        {        
//            BACOM(startSample, bacomPara.getNormalSamples().get(startSample), bacomPara.getTumorSamples().get(startSample), bacomPara.getInputCdfFile(),
//                    bacomPara.getAnnotationFile(), bacomPara.getResultFolder(),  bacomPara.getSampleIds().get(startSample),  bacomPara.getGeneInterested());
//        }
//        else
//        {
//            int midSample = (startSample+endSample)/2;
//            
//            invokeAll(new CorrectNormalContam(bacomPara, startSample, midSample),
//                  new CorrectNormalContam(bacomPara,  midSample+1, endSample));
//            
//        }
        /*
         * if serial computing is wanted, these codes should be commented out  
         */
        
    }
    
    
    public void BACOM(int loci, String normalSample, String tumorSample, String cdf, String annotation, String outputFileName, String sampleID, Gene geneInterested)
    {
        AffyCopyNumberData cnv;
	CopyNumberDetection cnvDetection;
	Noticcor D;
        
        Logger.logging("Start analyzing sample: "  + tumorSample + "..." );
        Logger.logging("Reading original copy nubmer signals from raw .CEL file...");
	cnv = new AffyCopyNumberData(normalSample, tumorSample, cdf, annotation);

	if (cnv.isSuccessfulOpenFile && cnv.isSuccessfulOpenAnnotationFile) {

            Logger.logging("Normalizing copy number signals...");
	    cnv.localNormalization();

            Logger.logging("Detecting copy number aberrations...");
	    cnvDetection = new CopyNumberDetection(cnv);

            Logger.logging("Estimating normal tissue contamination...");
	    D = new Noticcor(cnv, cnvDetection);

	    D.correctNTC();
            
            if (Double.isNaN(D.alpha))        //if D.alpha is NaN, set is as 0, and still use it.
                D.alpha = 0.0;
 
            bacomPara.getAlpha()[loci] = (double)Math.round(D.alpha*100000)/100000;          
                               
            // Determine deletion type for the interested gene //
            
            
            if (geneInterested != null){
                
                 Logger.logging("Detect the deletion type of interested genome region");
                 boolean isHomoDeletion = false;        
                 boolean isHemiDeletion = false;

                 for (int i = 0; i < D.deletionSegChrID.size(); i++) {
    //                   
                       String a = D.deletionSegChrID.get(i);
                       String b = geneInterested.getChr();
                       boolean c = (a.equals(b));

                       if ((D.deletionSegChrID.get(i)).equals(geneInterested.getChr()))
                       {
                           if (cnv.location[D.deletionSegStartLocation.get(i)] < geneInterested.getEnd())
                           {
                               if (cnv.location[D.deletionSegEndLocation.get(i)] > geneInterested.getStart())
                               {
                                   if (D.deletionSegProbHemiDeletion.get(i)< 0.01)
                                   {
                                       if (isHomoDeletion)
                                           continue;
                                       else{
                                            isHomoDeletion = true;
                                            GeneDeletionInfo.getInstance().addHomoDeletion();
                                       }
                                   }
                                   else if (D.deletionSegProbHomoDeletion.get(i) < 0.01)
                                   {
                                       if(isHemiDeletion)
                                           continue;
                                       else{
                                            isHemiDeletion = true;
                                            GeneDeletionInfo.getInstance().addHemiDeletion();
                                       }
                                   }
                               }
                           }
                       }
                   }
                }
               
            Logger.logging("Writing BACOM results to files...");

	    try {
                
                outputFileName = outputFileName + File.separator + "BACOMResults";
		if (!(new File(outputFileName)).exists()) {
		    boolean isSucessful = (new File(outputFileName)).mkdir();
		}                
                outputFileName = outputFileName + File.separator + sampleID + "_output";                
                writeResults(cnv, cnvDetection, D, outputFileName);
	    } catch (Exception e) {
                Logger.logging("Failed to create the file!" );
	    }

	} else {
	    if (!cnv.isSuccessfulOpenFile) {
		Logger.logging("\n  Errors in reading CEL/CDF files." +
				   "\n  Please check input files.\n");
	    } else if (!cnv.isSuccessfulOpenAnnotationFile) {
		Logger.logging("\n  Error in reading annotation files."
				   + "\n  Please make sure correct annotation files "
				   + "are placed under folder annotation.\n");
	    }
	}        
        cnv = null;        
        cnvDetection = null;        
        D = null;
    }
    
    
   private  void writeResults(AffyCopyNumberData cnv, CopyNumberDetection cnvDetection,
				     Noticcor D, String FileName) throws IOException {
               
	BufferedWriter resultFile = null;

	double[] segMean = new double[cnv.probeSetID.length];

        int k = 0;
        for (int i = 0; i < cnvDetection.chromosome.size(); i++) {

            for (int j = cnvDetection.segStartPos.get(i); j <= cnvDetection.segEndPos.get(i); j++) {
                segMean[k] = cnvDetection.segStatus.get(i);
                k ++;
            }
        }
                
        try {            

            resultFile = new BufferedWriter(new FileWriter(FileName + "data.csv"), 30*1024*1024);
            
            double correctCNA=0.0;
            
	    for (int i = 0; i < cnv.probeSetID.length; i ++) {
//		String entry = cnv.probeSetID[i] + "," + cnv.chrID[i] + "," + cnv.location[i] + ","
//		    + cnv.copyNumber[i] + "," + segMean[i];
                
                correctCNA = (segMean[i]-2*D.alpha)/(1-D.alpha);
                segMean[i] = (double) Math.round(segMean[i]*10000)/10000;
                correctCNA = (double) Math.round(correctCNA*10000)/10000;
                 cnv.copyNumber[i] = (double) Math.round( cnv.copyNumber[i]*10000)/10000;
                
                if (correctCNA<0)
                    correctCNA = Double.MIN_NORMAL;

                String entry = cnv.chrID[i] + "," + cnv.probeSetID[i] +"," + cnv.location[i] + ","
		    + cnv.copyNumber[i] + "," + segMean[i] + "," + correctCNA;
                
                resultFile.write(entry.toCharArray());
                
		resultFile.newLine();
            }
	}
        catch (Exception e){
            Logger.logging("Error writing BACOM results to files..." );
        }
        finally {
	    if (resultFile != null) {
                resultFile.close();
            }
            
        }

	BufferedWriter detectionResultFile = null;

        try {
            detectionResultFile = new BufferedWriter(new FileWriter(FileName + "detectionResult.csv"), 30*1024*1024);

	    for (int i = 0; i < cnvDetection.chromosome.size(); i ++) 
            {
		String entry = cnvDetection.chromosome.get(i) + ","
		    + cnvDetection.segStartPos.get(i) + ","
		    + cnvDetection.segEndPos.get(i) + ","
		    + cnvDetection.segStatus.get(i);
		detectionResultFile.write(entry.toCharArray());
		detectionResultFile.newLine();
            }
            
	} finally {
	    if (detectionResultFile != null) {
                detectionResultFile.close();
            }
        }


	BufferedWriter correctionResultFile = null;

        try {
            correctionResultFile = new BufferedWriter(new FileWriter(FileName + "BACOMresults.txt"), 30*1024*1024);
            String entry = D.alpha + "";
	    correctionResultFile.write(entry.toCharArray());
	    correctionResultFile.newLine();

	} finally {
	    if (correctionResultFile != null) {
                correctionResultFile.close();
            }
        }


        /*Write BACOM results with physical locations on the chromosomes   */

	correctionResultFile = null;

        try {
            correctionResultFile = new BufferedWriter(new FileWriter(FileName + "BACOM_locations.txt"), 30*1024*1024);
            
            String colName = "Chr" + "\t" + "from" + "\t" + "to" + "\t" + "Pr of hemi-deletion" + "\t" + "Pr of homo-deletion" + "\t" + "Normal Tissue fraction" ;
                        
            correctionResultFile.write(colName.toCharArray());
            correctionResultFile.newLine();

	    for (int i = 0; i < D.deletionSegChrID.size(); i ++) {

		String entry = D.deletionSegChrID.get(i) + "\t"
		    + cnv.location[D.deletionSegStartLocation.get(i)]+ "\t"
                    + cnv.location[D.deletionSegEndLocation.get(i)] + "\t"
		    + D.deletionSegProbHemiDeletion.get(i) + "\t"
		    + D.deletionSegProbHomoDeletion.get(i) + "\t"
		    + D.deletionSegNormalTissueFraction.get(i);

		correctionResultFile.write(entry.toCharArray());
		correctionResultFile.newLine();

            }
	} finally {
	    if (correctionResultFile != null) {
                correctionResultFile.close();
            }

       }

    }






}