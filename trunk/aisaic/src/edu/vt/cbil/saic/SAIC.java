/*
 * Package Name: SAIC "Significant Aberrations in Cancer Genomes"
 * Function: To detect the copy number alterations that are associated with cancer.
 * Author: Xuchu(Bella) HOU
 * Time:  Sep. 2012
 */
package edu.vt.cbil.saic;

import edu.vt.cbil.util.Logger;
import edu.vt.cbil.util.GeneSet;
import edu.vt.cbil.util.SaicParameters;
import java.lang.*;
import java.util.*;
import java.io.*;
import java.util.concurrent.RecursiveAction;  //Concurrent computing
import java.util.concurrent.ForkJoinPool;


/**
 * @code    SAIC detect Significant Copy number Aberrations (SCAs) 
 * in cancer genomes. 
 * @author  Xuchu(Bella) HOU
 * @version 2012.09
 */

public class SAIC extends RecursiveAction{

    
     /**
     * input parameters for SAIC analysis
     */
    
    
    private SaicParameters SAICPara;      
    
    /**
     * copy number alteration type, could be +1, 0 and -1 
     * +1 means amplification; 0 means needing concurrent computing
     * -1 means deletion
     */
    
    private int AlterType;    
 
    /**
     * minimum p-value when there is no permutated Uscore larger than observed Uscore
     */
    
    private  double MinPvalue;       
    
    /**
     * number of probes in the cancer genome
     */

    private  int inputProbeSize;    
    
    /**
     * number of samples
     */
    
    private  int inputSampleSize;  
    
    /**
     * number of Copy Number Alteration Probes (CNAP)
     */
    private  int pSize;            
    
    /**
     * number of Copy Number Alteration Units/Regions (RSize)
     */
    
    private  int RSize;             
    
    /**
     * number of permutation times
     */
    
    private  int PermuteSize;  
    
     /**
     * headers of the matrix, including chromosome, probe names, and location info
     */
    
    public  String headers[][]; 
    
    /**
     * input data matrix, each row corresponds to a single probe
     * each column corresponds to the copy number profile of a sample
     */
    
    private  double input_data[][];  
        
    
    /**
     * input amplification data matrix
     */
    
    private  double input_amp[][];    
    
    /**
     * input amplification data matrix
     */
    
    private  double input_del[][];    
    
    /**
     * p-values of the constructed CNA units
     */
    
    private  double pValue[];         
    
    /**
     * a ArrayList to store copy number alteration probes
     */    

    private ArrayList<Integer> CNAP;       
    
    /**
     * a ArrayList to store copy number alteration intervals
     */
    private ArrayList<locus> CNAI;         
    
    /**
     * a ArrayList to store copy number alteration units/regions
     */
    
    private ArrayList<locus> CNAR;        
    
    /**
     * a backup for copy number alteration units/regions 
     */
    
    private ArrayList<locus> CNAR_Copy;
    
    /**
     * a ArrayList to store Uscores for CNA units
     */
    
    private ArrayList<Double> UScore ;   
    
    /**
     * a backup for Uscores
     */
    
    private ArrayList<Double> UScore_Copy;
    
    /**
     * a ArrayList to store the lengths of the CNA units
     */
    
    private ArrayList<Integer> lenR ;   
    
    /**
     * a backup for lenR
     */
    
    private ArrayList<Integer> lenR_Copy;
    
    /**
     * a ArrayList to store the unique lengths of CNA units
     */
    
    private ArrayList<Integer> ulenR;
    
    /**
     * minimum size of CNA units
     */
               
    private static int MINCNARLen = 6;  
    
    /**
     * down-sampling rate
     */
    
    private static final int DownSamplingRate = 4;  
    

    /**
     * A constructor for the class SAIC.
     * 
     * @param para   Parameters for SAIC analysis
     * @param Type   Analysis type
     * @param input  input data matrix
     */
    
    public SAIC(SaicParameters para, int Type, double input[][] )
    {
        SAICPara = para;
        
        MinPvalue = 1.0/(1.0+SAICPara.getNp());
        
        AlterType = Type;
        
        input_data = input;
        
        if(SAICPara.isQuickSAIC() == true)
            
            MINCNARLen = MINCNARLen/DownSamplingRate;
    }
    
     /**
     * A constructor for the class SAIC, default concurrent computing
     * 
     * @param para   Parameters for SAIC analysis
     * @param input  input data matrix
     */

    public SAIC(SaicParameters para, double input[][])
    {
        SAICPara = new SaicParameters();
        
        SAICPara = para;
        
        MinPvalue = 1.0/(1.0+SAICPara.getNp());
        
        AlterType = SaicParameters.Concur;
        
        input_data = input;
        
        if(SAICPara.isQuickSAIC() == true)
            
            MINCNARLen = MINCNARLen/DownSamplingRate;
    }

    
    
    @Override
    protected void compute()
    {
        
          // Sequential analysis of Amplification and Deletion   
//          SCAsDetection(Parameters.AMP);
//          SCAsDetection(Parameters.DEL);
        
        // Parallel analysis of Amplification and Deletion
        
        if (AlterType == SaicParameters.Concur)
        {
            SaicParameters SAICParaAmp = new SaicParameters();
            SAICParaAmp.Copy(SAICPara);
            SaicParameters SAICParaDel = new SaicParameters();
            SAICParaDel.Copy(SAICPara);
            Logger.logging("Detecting amplified and deleted SCAs parallelly...");
            System.out.println("Detecting amplified and deleted SCAs parallelly...");
            invokeAll(new SAIC(SAICParaAmp, SaicParameters.AMP, input_data),
                    new SAIC (SAICParaDel, SaicParameters.DEL, input_data));
        }
        else
        {                    
            SCAsDetection(AlterType);
        }

    }

    /**
     * SCAs detection
     * @param Type, could be amplification and deletion 
     */
    
    public void SCAsDetection (int Type)
    {
        if (Type == SaicParameters.AMP)
        {       
            Logger.logging("Detecting amplifed copy number probes...");
            System.out.println("Detecting amplifed copy number probes...");
            GetCNAs(Type);            
            CNAP_Test(input_amp);    
            Logger.logging("There are totally " + pSize + "amplified copy number probes.");
            System.out.println("There are totally " + pSize + "amplified copy number probes.");
            if (pSize !=0)
            {                
                Logger.logging("Constructing amplifed copy number units...");
                CNAP_merge();
                CNAI_divide(input_amp);
                Logger.logging("There are totally " + RSize + " amplified copy number units");
                Logger.logging("Test the significance of amplifed copy number units...");
                System.out.println("There are totally " + RSize + " amplified copy number units");
                System.out.println("Test the significance of amplifed copy number units...");
                CNAR_Test(input_amp);

                /*  Free unused memory */
                CNAR = null;                
                UScore = null;
                lenR = null;
                ulenR = null;
                Logger.logging("Writing out the results...");
                System.out.println("Writing out the results...");
                WriteFile(Type);          //Here AMP=true, DEL = False
                input_amp = null;
            }
            else{
                Logger.logging("There is no copy number amplification probes, this is probably due to too few samples.");            
                System.out.println("There is no copy number amplification probes, this is probably due to too few samples.");
            }
            
            CNAR_Copy = null;
            lenR_Copy = null;
            UScore_Copy = null;
        }
        else
        {
            Logger.logging("Detecting deleted copy number probes...");
            System.out.println("Detecting deleted copy number probes...");
            GetCNAs(Type);
            CNAP_Test(input_del);
            Logger.logging("There are totally " + pSize + " deleted copy number probes.");
            System.out.println("There are totally " + pSize + " deleted copy number probes.");
            if (pSize !=0)
            {
                Logger.logging("Constructing deleted copy number units...");
                System.out.println("Constructing deleted copy number units...");
                CNAP_merge();
                CNAI_divide(input_del);
                Logger.logging("There are totally " + RSize + " deleted copy number units");
                System.out.println("There are totally " + RSize + " deleted copy number units");
                Logger.logging("Test the significance of deleted copy number units..."); 
                System.out.println("Test the significance of deleted copy number units..."); 
                CNAR_Test(input_del);

                /* Free unused memory */
                CNAR = null;                
                UScore = null;
                lenR = null;
                ulenR = null;

                Logger.logging("Writing out the results for deleted SCA detection...");
                System.out.println("Writing out the results for deleted SCA detection...");
                WriteFile(Type);          //Here AMP=true, DEL = False
                
                input_del = null;
                CNAR_Copy = null;
                lenR_Copy = null;
                UScore_Copy = null;
            }
            else{
                Logger.logging("There is no copy number deletion probes, this is probably due to too few samples.");
                System.out.println("There is no copy number deletion probes, this is probably due to too few samples.");
            }
        }
    }

    /**
     * get copy number amplification/deletion matrix
     * @param Type   indicating the copy number alteration type
     */
  public void GetCNAs(int Type)
    {
        
        int i, j;

        double temp;

        if (SAICPara.isQuickSAIC() == false)
        {
            inputProbeSize = input_data.length;
            inputSampleSize = input_data[0].length;
            headers = new String [inputProbeSize][3];

        if (Type == SaicParameters.AMP)
        {
            input_amp = new double[inputProbeSize][inputSampleSize];

            for (i=0; i<inputProbeSize; i++)
            {
                headers[i][0] = DetectSCAs.chr[i];
                headers[i][1] = DetectSCAs.markerName[i];
                headers[i][2] = DetectSCAs.markerID[i];                

                for (j=0; j<inputSampleSize; j++)
                {
                    temp = input_data[i][j];
                    if (temp > SAICPara.getTa())
                        input_amp[i][j] = temp;
                    else
                        input_amp[i][j] = 0;
                }
            }
        }
        else
        {
            input_del = new double[inputProbeSize][inputSampleSize];

            for (i=0; i<inputProbeSize; i++)
            {                
                headers[i][0] = DetectSCAs.chr[i];
                headers[i][1] = DetectSCAs.markerName[i];
                headers[i][2] = DetectSCAs.markerID[i]; 

                for (j=0; j<inputSampleSize; j++)
                {
                    temp = input_data[i][j];
                    
                    if (temp < SAICPara.getTd())
                        input_del[i][j] = temp;
                    else
                        input_del[i][j] = 0;
                }            
            }
        }
        }
        else
        {
            inputProbeSize = input_data.length/DownSamplingRate;
            inputSampleSize = input_data[0].length;
            headers = new String [inputProbeSize][3];
            if (Type == SaicParameters.AMP)
            {
            input_amp = new double[inputProbeSize][inputSampleSize];

            for (i=0; i<inputProbeSize; i++)
            {
                
                headers[i][0] = DetectSCAs.chr[DownSamplingRate*i];
                headers[i][1] = DetectSCAs.markerName[DownSamplingRate*i];
                headers[i][2] = DetectSCAs.markerID[DownSamplingRate*i]; 

                for (j=0; j<inputSampleSize; j++)
                {
                    temp = input_data[DownSamplingRate*i][j];
                    if (temp > SAICPara.getTa())
                        input_amp[i][j] = temp;
                    else
                        input_amp[i][j] = 0;
                }
            }
            }
        else
        {
            input_del = new double[inputProbeSize][inputSampleSize];

            for (i=0; i<inputProbeSize; i++)
            {
                
                headers[i][0] = DetectSCAs.chr[DownSamplingRate*i];
                headers[i][1] = DetectSCAs.markerName[DownSamplingRate*i];
                headers[i][2] = DetectSCAs.markerID[DownSamplingRate*i];

                for (j=0; j<inputSampleSize; j++)
                {
                    temp = input_data[DownSamplingRate*i][j];

                    if (temp < SAICPara.getTd())
                        input_del[i][j] = temp;
                    else
                        input_del[i][j] = 0;
                }
            }
        }

        }
    }
  
  /**
   * detect copy number alteration probes.
   * @param CNAMatrix   copy number alteration matrix
   */

    public void CNAP_Test(double CNAMatrix[][])  // Call the CNA probes
    {
        CNAP = new ArrayList<Integer>(inputProbeSize);
        for (int i=0; i<inputProbeSize; i++)
        {
            for (int j=0; j<inputSampleSize; j++)
            {
                if (CNAMatrix[i][j]!=0)
                {
                    CNAP.add(i);
                    break;
                }
            }
        }
        pSize = CNAP.size();   //Can't use CNAP.capacity, which is the initial capacity, .size() is the number of elements
        PermuteSize = pSize;
        CNAP.trimToSize();      // Trims the capacity of this vector to be the vector's current size

    }
    
    /**
     * group the consecutive CNA probes into one interval
     */

    public void CNAP_merge()       //merge CNA probes into CNA intervals
    {
      //  CNAI = new ArrayList<locus>(200);  //Set an initial size of 200
        CNAI = new ArrayList<locus>(pSize/100);

        locus temp = new locus();
        int i, index;       

       temp.start = CNAP.get(0).intValue();  //CNAP need to be declared as Integer vector,
                                              //otherwise the return type of CNAP.elementAt(i) would be Object
       temp.end = CNAP.get(0).intValue();

       if (pSize>1)
       {
            for (i=1; i<pSize; i++)
            {
                index = (CNAP.get(i-1).intValue())+1;
                if (CNAP.get(i)!=index)
                {
                    temp.end = CNAP.get(i-1).intValue();
                    if ((temp.end - temp.start +1) > MINCNARLen)
                        CNAI.add(temp);
                    temp = new locus();
                    temp.start = CNAP.get(i);
                }              
            }
            temp.end = CNAP.get(pSize-1).intValue();
       }
       CNAI.add(temp);
       CNAI.trimToSize();  // Trims the capacity of this vector to be the vector's current size
       CNAP = null;        // Release the memory
    }

    
    /**
     * group the highly correlated CNA probes in the CNA intervals into CNA units
     * @param CNAMatrix  copy number alteration matrix
     */

    public void CNAI_divide(double CNAMatrix[][])     // Get the CNA units
    {
    //  CNAR = new ArrayList<locus>(300);
    //  lenR = new ArrayList<Integer>(300);

      CNAR = new ArrayList<locus>(pSize/50);
      lenR = new ArrayList<Integer>(pSize/50);
      
      double CorrCoef;
      int i,j, temp_start, temp_end, len=0;
      int CNAI_num = CNAI.size();  //size of CNAI must be larger than 0
      locus temp = new locus();
      
      for (i=0; i<CNAI_num; i++)
      {
        temp_start = CNAI.get(i).start;
        temp_end = CNAI.get(i).end;
        temp.start = temp_start;
        temp.end = temp_end;
        if (temp_end > temp_start)
        {
            for (j = temp_start; j<temp_end; j++)
            {
                CorrCoef = COEF(j, j+1, CNAMatrix);
                if (CorrCoef < SAICPara.getTci())        // If CorrCoef<Th_Coef, then a break point appears
                {
                    temp.end = j;
                    len = temp.end - temp.start+1;
                    if (len>MINCNARLen)          //Won't store CNA units with length smaller than MINCNARLen
                    {
                        lenR.add(len);
                        CNAR.add(temp);
                    }
                    temp = new locus();
                    temp.start = j+1;
                }
            }
            temp.end = j;
        }
        len = temp.end - temp.start+1;
        if (len>MINCNARLen)
        {
            CNAR.add(temp);
            lenR.add(len);  // Score the lengths of CNA units
        }
        temp = new locus();
      }
      RSize = CNAR.size();
      CNAR.trimToSize();
      lenR.trimToSize();
      CNAI = null;        //Free the memory, CNAI won't be used in the subsequent operations
    }
    
    /**
     * calculate the Pearson correlation between consecutive probes
     * @param probe1      first probe 
     * @param probe2      second probe
     * @param CNAMatrix   copy number alteration matrix
     * @return            Pearson Correlation Coefficient between two CNA prbes  
     */

    public double COEF(int probe1, int probe2, double CNAMatrix[][])    //Calculate the correlation between consecutive probes
    {
        double coef, meanP1, meanP2, stdP1, stdP2, cov;
        coef=0;  
        int i;

        meanP1 =0; meanP2=0;
        for (i=0; i<inputSampleSize; i++)
        {
           meanP1 = meanP1+CNAMatrix[probe1][i];
           meanP2 = meanP2+CNAMatrix[probe2][i];
        }
        meanP1 = meanP1/(inputSampleSize);
        meanP2 = meanP2/(inputSampleSize);
        
        stdP1=0; stdP2=0; cov=0;
        for (i=0; i<inputSampleSize; i++)
        {
            stdP1 = stdP1+(CNAMatrix[probe1][i]-meanP1)*(CNAMatrix[probe1][i]-meanP1);
            stdP2 = stdP2+(CNAMatrix[probe2][i]-meanP2)*(CNAMatrix[probe2][i]-meanP2);
            cov = cov+(CNAMatrix[probe1][i]-meanP1)*(CNAMatrix[probe2][i]-meanP2);
        }
        stdP1 = stdP1/(inputSampleSize-1);
        stdP2 = stdP2/(inputSampleSize-1);

        if (stdP1 ==0 ||stdP2==0)
            return (0);

        coef = cov/(Math.sqrt(stdP1)*Math.sqrt(stdP2)*(inputSampleSize-1));

        return coef;
    }
    
    /**
     * calculate the UScore of observed CNA units
     * @param CNAMatrix  copy number alteration matrix
     */

    public void CalUScore (double CNAMatrix[][])    //Calculate the U-score
    {
        int i,j,k,m;
        int start, end, len;
        double UScore_Temp;

        UScore = new ArrayList<Double>(RSize);

        for (i=0; i<RSize; i++)
        {
            start = CNAR.get(i).start;
            end = CNAR.get(i).end;
            len = lenR.get(i);
            UScore_Temp = 0;           

            for (j=0; j<inputSampleSize; j++)
            {
                for(k=start; k<=end; k++)
                {
                    UScore_Temp = UScore_Temp+CNAMatrix[k][j];
                 }
            }
            //  No need to compute Rscore for amp and del separately, just take the abs()
            UScore_Temp = Math.abs(UScore_Temp/(inputSampleSize*len));
            UScore.add(UScore_Temp);
        }
       // PermuteSize = sum_len;

    }

    /**
     * Find the max Uscore in the permutated matrix for each CNA unit length
     * @param PermuteUScore     Uscore for permuted CNA units
     * @param PermutedCNAs      permuted CNA units
     */

    public void FindMaxUScore (double PermuteUScore [], double PermutedCNAs[][])  //Find the maximum of permutated UScore
    {
        double MaxUScore;
        int i, j,k, Endpos;
        double UScore_temp=0.0;
        double ProbeScore[] = new double[PermuteSize];
        
        for (i=0; i< PermuteSize; i++)   // Calculate the score for each probeset
        {
            UScore_temp=0;
            for (j=0; j<inputSampleSize; j++)
            {
                UScore_temp = UScore_temp + PermutedCNAs[i][j];
            }
            ProbeScore[i] = Math.abs(UScore_temp/inputSampleSize);            
        }

        int len = ulenR.get(0);  // The first unit length is calculated
        Endpos = PermuteSize-len;
        double RegionMean[] = new double[Endpos];
        MaxUScore = 0;
        for (j=0; j<Endpos; j++)
        {
            UScore_temp=0;

            for (k=j; k<j+len; k++)
            {
                UScore_temp = UScore_temp + ProbeScore[k];
            }
            UScore_temp = Math.abs(UScore_temp/len);         //Calculate the UScore for permuted data
            RegionMean[j] = UScore_temp;

            if(UScore_temp>MaxUScore)
            {
                MaxUScore = UScore_temp;       //Find the maximum UScore in the permutated data
            }
         }
        PermuteUScore[0] = MaxUScore;

 // For units other than the first one, we could use the result of former units to reduce the calculation load.

        for (i=1; i<ulenR.size(); i++)
        {
            len = ulenR.get(i);
            Endpos = PermuteSize-len-1;
            MaxUScore = 0;

            for(j=0; j<Endpos; j++)
            {
                UScore_temp=0;
                for (k=j+ulenR.get(i-1); k<j+len; k++)
                {
                    UScore_temp = UScore_temp + ProbeScore[k];
                }
                UScore_temp = Math.abs(UScore_temp) + RegionMean[j]*ulenR.get(i-1);
                UScore_temp = UScore_temp/len;
                RegionMean[j] = UScore_temp;
                
                if(UScore_temp>MaxUScore)
                {
                   MaxUScore = UScore_temp;       //Find the maximum UScore in the permutated data
                }
             }   
            PermuteUScore[i] = MaxUScore;
        }
    }
    
    /**
     * exclude detected SCAs until there is no SCAs detected 
     * @param    PermuteUScore
     * @return   whether to continue the next iteration
     */

    public boolean SCAExclude(double PermuteUScore[][])   //SCA-excluding
    {
        double SigValue;
        int i,j, len, index;
        int NewPermuteSize;
        boolean flag = false;
        int MaxLen = Collections.max(ulenR);  // To make sure there are enough probes to consstruct
                                                 //the CNA unit with the largest length

        for (i=0; i<RSize; i++)
        {
            SigValue=0;
            index = ulenR.indexOf(lenR.get(i));

            for (j=0; j<SAICPara.getNp(); j++)
            {
                if(PermuteUScore[j][index]>=UScore.get(i))
                    SigValue++;
            }
            SigValue = Math.max(SigValue/SAICPara.getNp(), MinPvalue);

            if (SigValue<SAICPara.getTp())            //If finding SCAs, get rid of that CNAR and its length, and change RSize;
            {
                flag=true;
                NewPermuteSize = PermuteSize-lenR.get(i);
                if (NewPermuteSize>MaxLen)
                    PermuteSize = NewPermuteSize;
                CNAR.remove(i);
                UScore.remove(i);
                lenR.remove(i);
                RSize = CNAR.size();
                
            }
        }
        return flag;
    }

    /**
     * calculation of significance level for the CNA units
     * @param PermuteUScore   Uscores for permuted CNA units
     */
    
    public void SigCal (double PermuteUScore[][])      //Calculate the significance of CNARs
    {          
        int i,j, index;
        double SigValue;
        RSize = UScore_Copy.size();
        pValue = new double[RSize];
        
        for (i=0; i<RSize; i++)
        {
            SigValue=0.0;
            index = ulenR.indexOf(lenR_Copy.get(i).intValue());

            for (j=0; j<SAICPara.getNp(); j++)
            {
                if(PermuteUScore[j][index] >= UScore_Copy.get(i).doubleValue())
                    SigValue++;
            }
            pValue[i] = Math.max(SigValue/SAICPara.getNp(), MinPvalue);  // Two int division in java will retrun int!
        }
    }

    /**
     * Significance testing of the CNA units
     * @param CNAMatrix  copy number alteration matrix
     */
    
    public void CNAR_Test(double CNAMatrix[][])
    {
        int i,j,k,m;
      
        double MaxUScore[][];
        double pValue[];
        boolean LoopFlag=true;        
         
        lenR_Copy = new ArrayList(lenR);        
          
        CNAR_Copy = new ArrayList(CNAR);  // Significant CNAR will be removed, so make a copy for record
        
        CalUScore(CNAMatrix);                //Assign a Uscore to each CNA unit.
        
        UScore_Copy = new ArrayList(UScore);
        
        Set Uniqlen = new HashSet(lenR);    //Store the unique lengths in the lenR
        
        ulenR = new ArrayList(Uniqlen);
        
        Collections.sort(ulenR);        //Increase sorting of the ulenR
        
        if (ulenR.size() >1 )
        {

            MaxUScore = new double[SAICPara.getNp()][ulenR.size()];

            Logger.logging("Estimating the null distribution using random Permutation...");
            System.out.println("Estimating the null distribution using random Permutation...");
            if (SAICPara.isIterFlag()==true)
            {
                Logger.logging("Using SCA-excluding permutation strategy to iteratively estimating the null distribution");
                System.out.println("Using SCA-excluding permutation strategy to iteratively estimating the null distribution");
                while (LoopFlag==true)
                {
                    Logger.logging("The permutation size is:" + PermuteSize);
                    System.out.println("The permutation size is:" + PermuteSize);
                    double PermuteMatrix[][] = new double[PermuteSize][inputSampleSize];

                    ForkJoinPool myPool = new ForkJoinPool();

                    for (i = 0; i<SAICPara.getNp(); i++)
                    {
                        Permutation myPermute = new Permutation(PermuteMatrix, CNAMatrix, 0, (inputSampleSize-1));
                        myPool.invoke(myPermute);
                        FindMaxUScore(MaxUScore[i], PermuteMatrix);
                        myPermute = null;
                    }
                    LoopFlag = SCAExclude(MaxUScore);
                } 
                SigCal(MaxUScore);
            }
            else
            {
                Logger.logging("SCA-excluding permutation strategy is not used");
                System.out.println("SCA-excluding permutation strategy is not used");
                double PermuteMatrix[][] = new double[PermuteSize][inputSampleSize];
                ForkJoinPool myPool = new ForkJoinPool();
                for (i = 0; i<SAICPara.getNp(); i++)
                {
                    Permutation myPermute = this.new Permutation(PermuteMatrix, CNAMatrix, 0, (inputSampleSize-1));
                    myPool.invoke(myPermute);
                    FindMaxUScore(MaxUScore[i], PermuteMatrix);
                    SigCal(MaxUScore);
                }
            }
        }
        else
        {
            System.out.println("Only one CNA unit!");
            pValue = new double[RSize]; 
            pValue[0] = 1;
        }
        
    /* If adopting down sampling, the lenR_Copy should recover the true number of probes included in the CNA units  */
                
        if (SAICPara.isQuickSAIC() == true)  
        {
            for (i = 0; i<lenR_Copy.size(); i++)
            {
                lenR_Copy.set(i, (DownSamplingRate * lenR_Copy.get(i)));
            }
        }

    }

    /**
     * write the results into files
     * @param flag  indicating amplification or deletion
     */
    
    public void geneQuery(List geneList)
    {  
        Logger.logging("Querying the genes covered by these detected SCAs...");
        System.out.println("Querying the genes covered by these detected SCAs...");
        int startSCA, endSCA, chrSCA, startChrIndex, endChrIndex;
        for (int i=0; i<RSize; i++)
         {
             if (pValue[i]<0.05)
             { 
                 chrSCA = Integer.parseInt(headers[CNAR_Copy.get(i).start][0]);
                 startSCA = Integer.parseInt(headers[CNAR_Copy.get(i).start][2]);
                 endSCA = Integer.parseInt(headers[CNAR_Copy.get(i).end][2]);

                 startChrIndex = GeneSet.getGeneChr().indexOf(chrSCA);
                 endChrIndex = GeneSet.getGeneChr().lastIndexOf(chrSCA);
                 
                 for (int j=startChrIndex; j<endChrIndex; j++)
                 {
                     if (GeneSet.getGeneEnd().get(j)> startSCA)
                     {
                         if(GeneSet.getGeneStart().get(j) < endSCA)
                            geneList.add(GeneSet.getGene().get(j));
                     }                    
                 }
             }
         }

    }   
    
    public void calDelPortion(double[] homoDeletionFlag, double[] hemiDeletionFlag)
    {
        Logger.logging("Estimating the proportion of deleted copy number units among all these samples...");
        System.out.println("Estimating the proportion of deleted copy number units among all these samples...");
        String tempFile;
        
        for (int num=0; num<SAICPara.getSampleID().size(); num++)
        {
            tempFile = SAICPara.getBACOMResultDir() + SAICPara.getSampleID().get(num) + "_outputBACOM_locations.txt";            
            String temp;
            String[] strSplit;            
            int delSegNum = 0;
            
            try{
                BufferedReader input = new BufferedReader (new FileReader (tempFile));
                temp = input.readLine();
                while(input.readLine()!= null)
                {
                    delSegNum++;
                }
            
            } catch(IOException e)
            {
                System.out.println("Can not find the BACOM deletion type results");
                Logger.logging("Can not find the BACOM deletion type results");
            }
            
            int[] chr = new int[delSegNum];
            int [] start = new int[delSegNum];
            int [] end = new int[delSegNum];
            double[] proHomo = new double[delSegNum];
            double[] proHemi = new double[delSegNum];
        
            try{
                
                BufferedReader input = new BufferedReader (new FileReader (tempFile));                
                temp = input.readLine();    //Skip the first row, which is the column names                
                int k=0;
                while((temp = input.readLine())!= null)
                {
                    strSplit = temp.split("\t");
                    chr[k] = Integer.parseInt(strSplit[0]);
                    start[k] = Integer.parseInt(strSplit[1]);
                    end[k] = Integer.parseInt(strSplit[2]);
                    proHomo[k] = Double.parseDouble(strSplit[3]);
                    proHemi[k] = Double.parseDouble(strSplit[4]);
                    k++;                    
                }
                
            }catch(IOException e)
            {
               System.out.println("Can not find the BACOM deletion type results");
                Logger.logging("Can not find the BACOM deletion type results");
            }            
        
        boolean isHemiDeletion, isHomoDeletion;

        for (int i=0; i<RSize; i++)
        {
            isHemiDeletion = false;
            isHomoDeletion = false;
            
            for (int j = 0; j < delSegNum; j++) {
//                                  
                if (chr[j]==Integer.parseInt(headers[CNAR_Copy.get(i).start][0]))
                {
                    if (start[j] < Integer.parseInt(headers[CNAR_Copy.get(i).end][2]))
                    {
                        if (end[j] > Integer.parseInt(headers[CNAR_Copy.get(i).end][2]))
                        {
                            if (proHemi[j]< 0.01)
                            {
                                if (isHemiDeletion)
                                    continue;
                                else{
                                     isHemiDeletion = true;
                                     hemiDeletionFlag[i]+=1;
                                   }
                             }
                             else if (proHomo[j] < 0.01)
                             {
                                   if(isHomoDeletion)
                                       continue;
                                   else
                                   {
                                        isHomoDeletion = true;
                                        homoDeletionFlag[i]+=1;;
                                   }
                             }
                           }
                       }
                   }
               }
        }
    }
    }

    
    public void WriteFile(int flag)
    {
        String AmpResult, DelResult;
        List ampGeneList, delGeneList;

        if (flag == SaicParameters.AMP)
        {
            // Write out the amplified genes  //
            
            ampGeneList = new ArrayList<String>(5000);
            geneQuery(ampGeneList);            
            String geneListDir = SAICPara.getOutfile() + "_ampGene";            
            try{
                
                BufferedWriter geneList = new BufferedWriter(new FileWriter(geneListDir));                
                geneList.write("Genes Covered in amplified SCAs");
                geneList.newLine();
                for (int i=0; i<ampGeneList.size(); i++)
                {
                    geneList.write(ampGeneList.get(i) + "\t");                    
                    geneList.newLine();
                }
                
                geneList.close();             
                              
            }catch(IOException e)
            {
                System.out.println("File Writing Error!");
                Logger.logging("File Writing Error!");
            }
             
            // prepare data for plotting average copy number profile across samples.
            double[] ampAvgCNA = new double[inputProbeSize];
            int[] ampScaIndicator = new int[inputProbeSize];
            for (int i = 0; i< inputProbeSize; i++){
                double sum = 0.0;
                for (int j=0; j<inputSampleSize; j++){
                    sum = sum + input_amp[i][j];
                }
                ampAvgCNA[i] = sum/inputSampleSize;
            } 
            int start = 0, end = 0;
            for (int i=0; i<RSize; i++){
                if(pValue[i]<SAICPara.getTp()){
                    start = CNAR_Copy.get(i).start;
                    end = CNAR_Copy.get(i).end;
                    Arrays.fill(ampScaIndicator, start, end, 1);
                }
            }
            try{
                String avgCnaFile = SAICPara.getOutfile() + "_ampAvgCNA";
                BufferedWriter bufferAvgCNA = new BufferedWriter(new FileWriter(avgCnaFile));
                for (int i=0; i<inputProbeSize; i++){
                    bufferAvgCNA.write("" + ampAvgCNA[i] + "\t" + ampScaIndicator[i]);
                    bufferAvgCNA.newLine();
                }
                bufferAvgCNA.close();
                
            }catch (Exception e){
                System.out.println("File Writing Error!");
                Logger.logging("File Writing Error!");
            }
                
            // write out the detected amplified SCAs //
            AmpResult = SAICPara.getOutfile()  + "_ampSCA";        
            try{            
                BufferedWriter SAICBufferOut = new BufferedWriter(new FileWriter(AmpResult));
                String entry = "Chromosome" + "\t" + "Start Location" + "\t"
                     + "End Location" + "\t" + "Length" + "\t" + "U-Score" + "\t" + "P-value";
                SAICBufferOut.write(entry.toCharArray());
                SAICBufferOut.newLine();
                for (int i=0; i<RSize; i++)
                {
                 if (pValue[i]<1)
                 {                
                    entry = headers[CNAR_Copy.get(i).start][0] + "\t" + headers[CNAR_Copy.get(i).start][2] + "\t"
                            + headers[CNAR_Copy.get(i).end][2] + "\t" + lenR_Copy.get(i) + "\t"
                            + UScore_Copy.get(i) + "\t" + pValue[i];
                    
                    SAICBufferOut.write(entry.toCharArray());                
                    SAICBufferOut.newLine();
                 }
                }
                SAICBufferOut.close();
            }
            catch (IOException e3)
            {
                System.out.println("File Writing Error!");
                Logger.logging("File Writing Error!");
            }
            catch (NullPointerException e4)
            {
                System.out.println("File Writing Error!");
                Logger.logging("File Writing Error!");

            }
        }
        
        else if (flag == SaicParameters.DEL)
        {
            
             double[] hemiDeletionNum = new double[RSize];
             double[] homoDeletionNum = new double[RSize];   

             if (SAICPara.getSampleID() != null)
             {
                    calDelPortion (homoDeletionNum, hemiDeletionNum);

                    for (int i = 0; i<RSize; i++)
                    {
                        hemiDeletionNum[i] = hemiDeletionNum[i]/inputSampleSize;
                        homoDeletionNum[i] = homoDeletionNum[i]/inputSampleSize;
                    }
             }
            
             // Write out the deleted genes  //
            delGeneList = new ArrayList<String>(5000);
            geneQuery(delGeneList);            
            String geneListDir = SAICPara.getOutfile() + "_delGene";            
            try{
                
                BufferedWriter geneList = new BufferedWriter(new FileWriter(geneListDir));
                geneList.write("Genes Covered in deleted SCAs");
                geneList.newLine();
                for (int i=0; i<delGeneList.size(); i++)
                {
                    geneList.write(delGeneList.get(i) + "\t");                    
                    geneList.newLine();
                }                
                geneList.close();                                   
            }catch(IOException e)
            {
               Logger.logging("Error writing results!");
               System.out.println("Error writing results!");
            }
            
            // prepare data for plotting average copy number profile across samples.
            double[] delAvgCNA = new double[inputProbeSize];
            int[] delScaIndicator = new int[inputProbeSize];
            for (int i = 0; i< inputProbeSize; i++){
                double sum = 0.0;
                for (int j=0; j<inputSampleSize; j++){
                    sum = sum + input_del[i][j];
                }
                delAvgCNA[i] = sum/inputSampleSize;
            } 
            int start = 0, end = 0;
            for (int i=0; i<RSize; i++){
                if(pValue[i]<SAICPara.getTp()){
                    start = CNAR_Copy.get(i).start;
                    end = CNAR_Copy.get(i).end;
                    Arrays.fill(delScaIndicator, start, end, 1);
                }
            }
            try{
                String avgCnaFile = SAICPara.getOutfile() + "_delAvgCNA";
                BufferedWriter bufferAvgCNA = new BufferedWriter(new FileWriter(avgCnaFile));
                for (int i=0; i<inputProbeSize; i++){
                    bufferAvgCNA.write("" + delAvgCNA[i] + "\t" + delScaIndicator[i]);
                    bufferAvgCNA.newLine();
                }
                bufferAvgCNA.close();
                
            }catch (Exception e){
                System.out.println("Error writing file!");
                Logger.logging("Error writing file!");
            }
            
            // Write out the SCA detecting resutls.
            DelResult = SAICPara.getOutfile() + "_delSCA";
            
            if (SAICPara.getSampleID() == null)
            {
            
                try{

                    BufferedWriter SAICBufferOut = new BufferedWriter(new FileWriter(DelResult));

                    String entry = "Chromosome" + "\t" + "Start Location" + "\t" + "End Location" + "\t" + "Length" + 
                            "\t" + "U-Score" + "\t" + "P-value" + "\t" ;

                    SAICBufferOut.write(entry.toCharArray());             
                    SAICBufferOut.newLine();             
                    for (int i=0; i<RSize; i++)
                    {
                        if (pValue[i]<1)
                        {
                            entry = headers[CNAR_Copy.get(i).start][0] + "\t" + headers[CNAR_Copy.get(i).start][2] + "\t"
                            + headers[CNAR_Copy.get(i).end][2] + "\t" + lenR_Copy.get(i) + "\t"
                            + UScore_Copy.get(i) + "\t" + pValue[i] + "\t" ;

                        SAICBufferOut.write(entry.toCharArray());                
                        SAICBufferOut.newLine();
                        }
                    }                
                    SAICBufferOut.close();
                }
                catch (IOException e5)
                {
                    Logger.logging("File Writing Error!");
                    System.out.println("File Writing Error!");
                }
                catch (NullPointerException e6)
                {
                    Logger.logging("File Writing Error!");
                    System.out.println("File Writing Error!");
                }
            }
            else
            {
                try{

                    BufferedWriter SAICBufferOut = new BufferedWriter(new FileWriter(DelResult));

                    String entry = "Chromosome" + "\t" + "Start Location" + "\t" + "End Location" + "\t" + "Length" + 
                            "\t" + "U-Score" + "\t" + "P-value" + "\t" + "Proportion of HomoDeletion" + "\t" + "Proportion of HemiDeletion";

                    SAICBufferOut.write(entry.toCharArray());

                    SAICBufferOut.newLine();

                    for (int i=0; i<RSize; i++)
                    {
                        if (pValue[i]<1)
                        {
                            entry = headers[CNAR_Copy.get(i).start][0] + "\t" + headers[CNAR_Copy.get(i).start][2] + "\t"
                            + headers[CNAR_Copy.get(i).end][2] + "\t" + lenR_Copy.get(i) + "\t"
                            + UScore_Copy.get(i) + "\t" + pValue[i] + "\t" + homoDeletionNum[i] + "\t" + hemiDeletionNum[i];

                        SAICBufferOut.write(entry.toCharArray());

                        SAICBufferOut.newLine();
                        }
                    }

                    SAICBufferOut.close();
                }
                catch (IOException e5)
                {
                    Logger.logging("File Writing Error!");
                    System.out.println("File Writing Error!");
                }
                catch (NullPointerException e6)
                {
                   Logger.logging("File Writing Error!");
                    System.out.println("File Writing Error!");
                }

                homoDeletionNum = null;
                hemiDeletionNum = null;
            }
        }

    }
 
/**
 * locus 
 */
private class locus{                 //Define locus, used in CNA units
    
    /**
     * start location of a copy number region
     */
    int start;
    
    /**
     * end location of a copy number region
     */
    
    int end;

    /**
     * constructor of locus
     */
    public locus()
    {
        start=0;        
        end=0;
    }
}

/**
 * Permutation    within sample permutation, inherit from RecursiveAction, used to 
 */
private class Permutation extends RecursiveAction{

    /**
     *  permuted copy number alteration matrix
     */
    
    double PermuteCNAs[][];
    
    /**
     * original copy number alteration matrix;
     */
    
    double CNAMatrix[][];
    
    /**
     * start sample for a permutation thread
     */

    int StartSample;
    
    /**
     * end sample for a permutation thread
     */
    
    int EndSample;
    
    /**
     * number of samples permuted in a thread
     */
    
    int RegionSize; 

    public Permutation(double PermutedCNAs[][], double CNAMatrix[][], int start, int end)
    {
        this.PermuteCNAs = PermutedCNAs;
        this.CNAMatrix = CNAMatrix;
        this.EndSample = end;
        this.StartSample = start;
        this.RegionSize = lenR.size();

    }

    public void Permute( )  //Sample-wise permutation
    {
        int i,j,k,m;
        int index, MaxIndex, startP, pos, len;

        Vector<Integer> TempSet = new Vector<Integer>(RegionSize);   //Store the index of the CNA regions
        Vector<Integer> TempSet_Copy = new Vector<Integer>(RegionSize);
        for (i=0; i<RegionSize; i++)         //
        {
            TempSet.add(i);
        }

        Random rn  = new Random();
        for (j=StartSample; j<=EndSample; j++)         //Start of sample permutaion
        {
                MaxIndex = RegionSize;
                TempSet_Copy.addAll(TempSet);
                startP=0;
                for (k=0; k<RegionSize; k++)
                {
                    index = rn.nextInt(MaxIndex);                   // Randomly draw a index of the region IDs
                    pos = TempSet_Copy.elementAt(index);           // Get the corresponding posth CNA region
                    len = lenR.get(pos);
                    for (m=startP; m<(startP+len); m++)
                    {
                        PermuteCNAs[m][j] = CNAMatrix[CNAR.get(pos).start+m-startP][j];
                    }
                    MaxIndex--;
                    TempSet_Copy.remove(index);         // Remove the already drawed index, without replacement
                    startP = startP+len;                // Increase the start position for inserting new CNARegion
                }
                                    // End of sample permutation
         }

/*******************************************************************************************/
//**         The following block is used for testing the permutation       **//
//         try{
//         String out = "test/tempout.txt";
//         BufferedWriter BufferOut = new BufferedWriter(new FileWriter(out));
//         for (m=0; m<PermuteSize; m++)
//         {
//             String entry = new String();
//             for (k=0; k<inputSampleSize; k++)
//             {
//                entry = entry + "\t" + PermutedCNAs[m][k];
//             }
//            BufferOut.write(entry.toCharArray());
//            BufferOut.newLine();
//         }
//         BufferOut.close();
//        }
//        catch (IOException e3)
//        {
//            System.out.println("File Writing Error!");
//            e3.printStackTrace();
//        }
//**         The above block is used for testing the permutation       **//
/******************************************************************************************************/

    }

    protected void compute()
    {
        int NumSample = EndSample - StartSample +1;
        if (NumSample < 20)
        {
            Permute();
        }
        else
        {
            int SubEnd = StartSample + NumSample/2;
            int SubStart = SubEnd+1;
            invokeAll(new Permutation(PermuteCNAs, CNAMatrix, StartSample, SubEnd),
                    new Permutation(PermuteCNAs, CNAMatrix, SubStart, EndSample));
            
        }
    }

}
 

}

