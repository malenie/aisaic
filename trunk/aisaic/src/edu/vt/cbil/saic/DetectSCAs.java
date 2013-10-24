/*
 * Class Name: DetectSCAs
 * Function: Call SAIC (For concurrent computing use)
 * Author: Xuchu (Bella) HOU
 * Email: houxc0511@gmail.com // bella@vt.edu
 * Time:  Sep. 2012
 */

package edu.vt.cbil.saic;

import edu.vt.cbil.exceptions.TooSmallSampleSizeException;
import edu.vt.cbil.util.SaicParameters;
import java.io.*;
import java.util.concurrent.*;
import java.util.Vector;
import edu.vt.cbil.util.Logger;

/*  To analyze two chromosomes at the same time   */

public class DetectSCAs extends RecursiveAction{
//    
     SaicParameters SCAPara;
     String[] infileSet;
     String[] outfileSet;
     
     public static String[] chr;
     public static String[] markerName;
     public static String[] markerID;
     
    
    public DetectSCAs(SaicParameters UsrPara, String[] inputfiles, String[] outputfiles)
    {
        SCAPara = UsrPara;
        infileSet = inputfiles;
        outfileSet = outputfiles;         
    
    }
    
    public void DetectingSCA() throws TooSmallSampleSizeException
    {   
        
        SCAPara.setInfile(infileSet[0]);
        SCAPara.setOutfile(outfileSet[0]);   
             
        double InputData[][] = ReadFile(SCAPara.getInfile());
        
        ForkJoinPool myForkJoinPool = new ForkJoinPool();
               
        SAIC mySAIC = new SAIC(SCAPara, SaicParameters.Concur, InputData);
        
        myForkJoinPool.invoke(mySAIC);  
        
        mySAIC = null;
        myForkJoinPool = null;
      
     }
    
    protected void compute()
    {
        
        if ((infileSet.length > 1) && (outfileSet.length > 1))
        {
            String[] infileSet1 = new String[1];
            String[] infileSet2 = new String[1];
            String[] outfileSet1 = new String[1];
            String[] outfileSet2 = new String[1];      
            infileSet1[0] = infileSet[0];
            infileSet2[0] = infileSet[1];
            outfileSet1[0] = outfileSet[0];
            outfileSet2[0] = outfileSet[1];
//            
           SaicParameters SCAPara1 =new SaicParameters();
           SCAPara1.Copy(SCAPara);
           SaicParameters SCAPara2 = new SaicParameters();
           SCAPara2.Copy(SCAPara);         
           
           
            invokeAll(new DetectSCAs(SCAPara1, infileSet1,outfileSet1),
                    new DetectSCAs(SCAPara2, infileSet2,outfileSet2));
        }        
        else {
            try{
                 DetectingSCA();
            }catch (Exception e){
                System.out.println(e);
                return;
            }
        }
           

    }

    
 /*  In order to realize concurrent computing, "ReadFile" is put outside class SAIC
         otherwise, either the input file need to be read twice, or unable to adopt concurrent computing  */
    
public double[][] ReadFile(String InFile) throws TooSmallSampleSizeException
    {
        int NumRow, NumCol;
               
        double DataIn[][];
        
        try{             //It's OK to specify more than one exception
            
            BufferedReader SAICBufferIN = new BufferedReader (new FileReader(InFile));

            try{
                String StrTemp;
                
                NumRow=0;

                StrTemp = SAICBufferIN.readLine();
                
                String StrSplit[]= StrTemp.split("\t");
                
                NumCol = StrSplit.length-3;     //Calculate the number of samples
                
                if ((NumCol-3)<=1)
                {
                    System.out.println("Too small sample size! Please add more samples");
                    throw new TooSmallSampleSizeException("Too small sample size!");
                }

                NumRow++;

                while ((StrTemp = SAICBufferIN.readLine())!= null)  //After combing with BACOM
                {
                    NumRow++;          //Calculate the number of probes
                }                              // Define the size of vector before using is a better choice than using add
                SAICBufferIN.close();  // This closes stream for both BufferedReader and InputStreamReader.

                int i=0;
                
                BufferedReader SAICBuffer_R = new BufferedReader (new FileReader(InFile));

                DataIn = new double[NumRow][NumCol];
//                
                chr = new String[NumRow];
                markerName = new String[NumRow];
                markerID = new String[NumRow];

                System.out.println("Reading the input file...");
                
                while ((StrTemp = SAICBuffer_R.readLine())!= null)
                {
                    StrSplit= StrTemp.split("\t");
                    
                    chr[i] = StrSplit[0];
                    markerName[i] = StrSplit[1];
                    markerID[i] = StrSplit[2];
//                    
//                    System.out.println(NumCol);
//
//                    
                    for(int j=0; j<NumCol; j++)
                    {
                        DataIn[i][j] = Double.parseDouble(StrSplit[j+3]);

                    }
                    i++;
                }                
                SAICBuffer_R.close();  // Only need to close the outer wrapper.
                
               return DataIn;
             }
            catch (IOException e0)
            {
                Logger.logging("Error reading file!" + "\n" + "\n");  
                System.exit(1);
                return null;
            }
        }
        catch(FileNotFoundException e1)
        {
            Logger.logging("File is not found! Please check it again." + "\n" + "\n");  
            System.exit(1);
            return null;
        }
        catch (SecurityException e2)
        {
            Logger.logging("File access is denied! Please check it again." + "\n" + "\n");
            System.exit(1);
            return null;
        }
    } 

}
