/*
 * Class Name: Parameters
 * Function: To store the parameters for SAIC
 * Author: Xuchu (Bella) HOU 
 * Time:  Sep. 2012
 */
package edu.vt.cbil.util;

import edu.vt.cbil.exceptions.FileDoesNotExistException;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Xuchu(Bella) HOU
 */

/**
 * Parameters class is used to store the parameters for SAIC 
 * @version   2012.09
 */

public class SaicParameters {
    
    /**
     * Analysis mode
     */

   private int analysisMode;
    /**
     * threshold for detecting copy number amplification 
     */
    
   private double ta;
   
   /**
     * threshold for detecting copy number deletion 
     */
   
   private double td;
   
   /**
     * number of permutation  
     */
   
   private int np;
   
   /**
     * p_value cutoff for significance test
     */
   
   private double tp;
   
   /**
     * cutoff for pearson correlation coefficient when constructing CNA units
     */
   
   private double tci;
   
   /**
     * SCA-excluding indication
     */
   
   private boolean iterFlag;
   
   /**
     * whether adopt quick SAIC
     */
   
   private boolean QuickSAIC;
   
   /**
     * input file or folder directory
     */
   
   private String infile;
   
   /**
     * output file or folder directory
     */
   
   private String outfile;
   
   private String infileDir;
   
   private String outfileDir;
   
   /**
    *  flag to indicate whether AISAIC is chosen 
    */
   
   private boolean isAISAIC;
   
   /**
    *  directory of the gene reference 
    */
   
   private String geneRefDir;
   
   /**
    *  the directory of BACOM results
    */
   
   private String BACOMResultDir;
   
   /**
    * sample ID of the data set
    */
   
   private ArrayList<String> sampleID;
   

   
   /**
     * static field amplification indicator
     */
   
   
   public static final int AMP = 1;
   
   /**
     * static field deletion indicator
     */
   
   public static final int DEL = -1;
   
   /**
     * static field concurrent computing indicator
     */
   
   public static final int Concur = 0;
   
   /**
    * Parameters class constructor
    */
   
    public SaicParameters()
    {
       ta = 0.322;      
       td = -0.415;       
       np = 1000;       
       tp = 0.0476;       
       tci = 0.9;       
       iterFlag=true;      
       QuickSAIC = true;      
       isAISAIC = false;            
    }

    /**
     * Parameters class constructor
     * @param in
     * @param out 
     */
    public SaicParameters(String in, String out)
    {
       ta = 0.322;      
       td = -0.415;      
       np = 1000;      
       tp = 0.0476;       
       tci = 0.9;       
       iterFlag=true;      
       infile = in;      
       outfile = out;       
       QuickSAIC = true;       
       isAISAIC = false;      
    }
    
    /**
     * Parameters class constructor
     * @param amp
     * @param del
     * @param numP
     * @param Pvalue_th
     * @param Corcoef_th
     * @param flag
     * @param qSAIC
     * @param in
     * @param out 
     */
    public SaicParameters(double amp, double del, int numP, double Pvalue_th, 
            double Corcoef_th, boolean flag, boolean qSAIC, String in, String out)
    {
        this.ta = amp;       
        this.td = del;       
        this.np = numP;        
        this.tp = Pvalue_th;        
        this.tci = Corcoef_th;        
        this.iterFlag = flag;        
        this.infile = in;        
        this.outfile = out;        
        this.QuickSAIC = qSAIC;        
        isAISAIC = false;      
        geneRefDir = null;      
        BACOMResultDir = null;      
        sampleID = null;      
    }
    
    /**
     * Parameters class copy
     * @param Origin 
     */
    
    public void Copy(SaicParameters Origin)
    {
        this.ta = Origin.ta;       
        this.td = Origin.td;        
        this.np = Origin.np;        
        this.tp = Origin.tp;        
        this.tci = Origin.tci;       
        this.iterFlag = Origin.iterFlag;        
        this.infile = Origin.infile;        
        this.outfile = Origin.outfile;        
        this.QuickSAIC = Origin.QuickSAIC;        
        this.geneRefDir = Origin.geneRefDir;        
        this.isAISAIC = Origin.isAISAIC;      
        this.geneRefDir = Origin.geneRefDir;      
        this.BACOMResultDir = Origin.BACOMResultDir;      
        this.sampleID = Origin.sampleID;
    }

    public void setTa(double ta) {
        this.ta = ta;
    }

    public void setTd(double td) {
        this.td = td;
    }

    public void setNp(int np) {
        this.np = np;
    }

    public void setTp(double tp) {
        this.tp = tp;
    }

    public void setTci(double tci) {
        this.tci = tci;
    }

    public void setIterFlag(boolean iterFlag) {
        this.iterFlag = iterFlag;
    }

    public void setQuickSAIC(boolean QuickSAIC) {
        this.QuickSAIC = QuickSAIC;
    }

    public void setInfile(String infile) {
        this.infile = infile;
    }

    public void setOutfile(String outfile) {
        this.outfile = outfile;
    }

    public void setIsAISAIC(boolean isAISAIC) {
        this.isAISAIC = isAISAIC;
    }

    public void setGeneRefDir(String geneRefDir) throws FileDoesNotExistException{
        File geneRefFile = new File(geneRefDir);
        if(!geneRefFile.exists())
            throw new FileDoesNotExistException("Can not find the specified file, please check again.");
        this.geneRefDir = geneRefDir;
    }

    public void setBACOMResultDir(String BACOMResultDir) {
        this.BACOMResultDir = BACOMResultDir;
    }

    public void setSampleID(ArrayList<String> sampleID) {
        this.sampleID = sampleID;
    }

    public double getTa() {
        return ta;
    }

    public double getTd() {
        return td;
    }

    public int getNp() {
        return np;
    }

    public double getTp() {
        return tp;
    }

    public double getTci() {
        return tci;
    }

    public boolean isIterFlag() {
        return iterFlag;
    }

    public boolean isQuickSAIC() {
        return QuickSAIC;
    }

    public String getInfile() {
        return infile;
    }

    public String getOutfile() {
        return outfile;
    }

    public boolean isIsAISAIC() {
        return isAISAIC;
    }

    public String getGeneRefDir() {
        return geneRefDir;
    }

    public String getBACOMResultDir() {
        return BACOMResultDir;
    }

    public ArrayList<String> getSampleID() {
        return sampleID;
    }

    public static int getAMP() {
        return AMP;
    }

    public static int getDEL() {
        return DEL;
    }

    public static int getConcur() {
        return Concur;
    }

    public String getInfileDir() {
        return infileDir;
    }

    public void setInfileDir(String infileDir) {
        this.infileDir = infileDir;
    }

    public String getOutfileDir() {
        return outfileDir;
    }

    public void setOutfileDir(String outfileDir) {
        this.outfileDir = outfileDir;
    }

    public int getAnalysisMode() {
        return analysisMode;
    }

    public void setAnalysisMode(int analysisMode) {
        this.analysisMode = analysisMode;
    }
    
    
    
 
}

    
