/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.vt.cbil.util;

import java.util.ArrayList;
import java.io.*;
import java.util.List;

class Chromosome{

    static final int LENCHR1 = 146401;
    static final int LENCHR2 = 153663;
    static final int LENCHR3 = 127766;
    static final int LENCHR4 = 120296;
    static final int LENCHR5 = 115672;
    static final int LENCHR6 = 112825;
    static final int LENCHR7 = 100996;
    static final int LENCHR8 = 98277;
    static final int LENCHR9 = 82168;
    static final int LENCHR10 = 93592;
    static final int LENCHR11 = 89525;
    static final int LENCHR12 = 87321;
    static final int LENCHR13 = 66067;
    static final int LENCHR14 = 57103;
    static final int LENCHR15 = 53556;
    static final int LENCHR16 = 54182;
    static final int LENCHR17 = 46632;
    static final int LENCHR18 = 52093;
    static final int LENCHR19 = 30299;
    static final int LENCHR20 = 43628;
    static final int LENCHR21 = 25111;
    static final int LENCHR22 = 24484;  
    static final int NUMPROBE = 1781657;

}
public class DataConversion {
    
    double[] alpha;    
    String bacomResultDir; 
    String saicInputDir;
    ArrayList<String> sampleID; 
    
    public DataConversion(double[] alpha, ArrayList<String> sampleID, String bacomResultDir, String saicInputDir ) 
    {
        this.alpha = alpha;
        this.sampleID = sampleID;  
        this.bacomResultDir = bacomResultDir;
        this.saicInputDir = saicInputDir;
    } 
    
    public void transform()throws FileNotFoundException, IOException
    {
        int numSample = alpha.length;
        
        int[] chrID = new int[Chromosome.NUMPROBE];
        
        String[] probeName = new String[Chromosome.NUMPROBE];
        
        int[] probeLoci = new int[Chromosome.NUMPROBE];
        
        double[][] chr1 = new double[Chromosome.LENCHR1][numSample];
        
        double[][] chr2 = new double[Chromosome.LENCHR2][numSample];
        
        double[][] chr3 = new double[Chromosome.LENCHR3][numSample];
        
        double[][] chr4 = new double[Chromosome.LENCHR4][numSample];
        
        double[][] chr5 = new double[Chromosome.LENCHR5][numSample];
        
        double[][] chr6 = new double[Chromosome.LENCHR6][numSample];
        
        double[][] chr7 = new double[Chromosome.LENCHR7][numSample];
        
        double[][] chr8 = new double[Chromosome.LENCHR8][numSample];
        
        double[][] chr9 = new double[Chromosome.LENCHR9][numSample];
        
        double[][] chr10 = new double[Chromosome.LENCHR10][numSample];
        
        double[][] chr11 = new double[Chromosome.LENCHR11][numSample];
        
        double[][] chr12 = new double[Chromosome.LENCHR12][numSample];
        
        double[][] chr13 = new double[Chromosome.LENCHR13][numSample];
        
        double[][] chr14 = new double[Chromosome.LENCHR14][numSample];
        
        double[][] chr15 = new double[Chromosome.LENCHR15][numSample];
        
        double[][] chr16 = new double[Chromosome.LENCHR16][numSample];
        
        double[][] chr17 = new double[Chromosome.LENCHR17][numSample];
        
        double[][] chr18 = new double[Chromosome.LENCHR18][numSample];
        
        double[][] chr19 = new double[Chromosome.LENCHR19][numSample];
        
        double[][] chr20 = new double[Chromosome.LENCHR20][numSample];
        
        double[][] chr21 = new double[Chromosome.LENCHR21][numSample];
        
        double[][] chr22 = new double[Chromosome.LENCHR22][numSample];
        
        double[][] wholeGenome = new double[Chromosome.NUMPROBE][numSample];
        

        String temp;
        
        String[] tempStringSplit;
        
        int col=0;        

        BufferedReader input = new BufferedReader(new FileReader(bacomResultDir+sampleID.get(0) + "_outputdata.csv"));     
        System.out.println("Reading copy number profile of " + bacomResultDir+sampleID.get(0));
        Logger.logging("Reading copy number profile of " + bacomResultDir+sampleID.get(0));

        int probeChr1=0, probeChr2=0, probeChr3=0, probeChr4=0, probeChr5=0, probeChr6=0, probeChr7=0,
                probeChr8=0, probeChr9=0, probeChr10=0, probeChr11=0, probeChr12=0, probeChr13=0, probeChr14=0,
                probeChr15=0, probeChr16=0, probeChr17=0, probeChr18=0, probeChr19=0, probeChr20=0,
                probeChr21=0,probeChr22=0;

        int i=0;

        while(i<Chromosome.NUMPROBE && (temp = input.readLine())!=null)
        {                                       
            tempStringSplit = temp.split(",");                       
            double tmp = Math.log(Double.parseDouble(tempStringSplit[5]))/Math.log(2) - 1;                    
            tempStringSplit[5] = "" + tmp;

            chrID[i] = Integer.parseInt(tempStringSplit[0]);
            probeName[i] = tempStringSplit[1];
            probeLoci[i] = Integer.parseInt(tempStringSplit[2]);                    
            wholeGenome[i][col] = Double.parseDouble(tempStringSplit[5]); 

            i++;

            switch(Integer.parseInt(tempStringSplit[0]))                        
            {
                case 1:
                    chr1[probeChr1][col] = Double.parseDouble(tempStringSplit[5]);                              
                    probeChr1 ++;
                    break;

                case 2:
                    chr2[probeChr2][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr2++;                                                  
                    break; 

                case 3:
                    chr3[probeChr3][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr3++;
                    break;                        

                case 4:
                    chr4[probeChr4][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr4++;
                    break; 

                case 5:
                    chr5[probeChr5][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr5++;
                    break; 

                case 6:
                    chr6[probeChr6][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr6++;
                    break; 

                case 7:
                    chr7[probeChr7][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr7++;
                    break; 

                 case 8:
                    chr8[probeChr8][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr8++;
                    break; 

                 case 9:
                    chr9[probeChr9][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr9++;
                    break; 

                 case 10:
                    chr10[probeChr10][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr10++;
                    break;

                 case 11:
                    chr11[probeChr11][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr11++;
                    break; 

                 case 12:
                    chr12[probeChr12][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr12++;
                    break; 

                 case 13:
                    chr13[probeChr13][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr13++;
                    break; 

                 case 14:
                    chr14[probeChr14][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr14++;
                    break; 

                 case 15:
                    chr15[probeChr15][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr15++;
                    break; 

                 case 16:
                    chr16[probeChr16][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr16++;
                    break; 

                 case 17:
                    chr17[probeChr17][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr17++;
                    break; 

                case 18:
                    chr18[probeChr18][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr18++;
                    break; 

               case 19:
                    chr19[probeChr19][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr19++;
                    break; 

               case 20:
                    chr20[probeChr20][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr20++;
                    break; 

               case 21:
                    chr21[probeChr21][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr21++;
                    break; 

              case 22:
                    chr22[probeChr22][col] = Double.parseDouble(tempStringSplit[5]);  
                    probeChr22++;
                    break; 
            }
        }
       
        col++;
        
        for (int iSample=1; iSample<numSample; iSample++)
        {            
            input = new BufferedReader(new FileReader(bacomResultDir+sampleID.get(iSample)+"_outputdata.csv"));
            System.out.println("Reading copy number profile of " + bacomResultDir+sampleID.get(iSample));
            Logger.logging("Reading copy number profile of " + bacomResultDir+sampleID.get(iSample));

            probeChr1=0; probeChr2=0; probeChr3=0; probeChr4=0; probeChr5=0; probeChr6=0; probeChr7=0;
                    probeChr8=0; probeChr9=0; probeChr10=0; probeChr11=0; probeChr12=0; probeChr13=0; probeChr14=0;
                    probeChr15=0; probeChr16=0; probeChr17=0; probeChr18=0; probeChr19=0; probeChr20=0;
                    probeChr21=0;probeChr22=0;

            i=0;
            while(i<Chromosome.NUMPROBE && (temp = input.readLine())!=null)
            {                    
                tempStringSplit = temp.split(",");                    
                double tmp = Math.log(Double.parseDouble(tempStringSplit[5]))/Math.log(2) - 1;                                        
                tempStringSplit[5] = "" + tmp;                    
                wholeGenome[i][col] = Double.parseDouble(tempStringSplit[5]);                    
                i++;                    
                switch(Integer.parseInt(tempStringSplit[0]))                        
                {
                    case 1:

                        chr1[probeChr1++][col] = Double.parseDouble(tempStringSplit[5]);  
                        break;


                    case 2:

                        chr2[probeChr2++][col] = Double.parseDouble(tempStringSplit[5]);  
                        break;                        


                    case 3:

                        chr3[probeChr3++][col] = Double.parseDouble(tempStringSplit[5]);  
                        break;     

                    case 4:

                        chr4[probeChr4++][col] = Double.parseDouble(tempStringSplit[5]);  
                        break;


                    case 5:

                        chr5[probeChr5++][col] = Double.parseDouble(tempStringSplit[5]);  
                        break;                        


                    case 6:

                        chr6[probeChr6++][col] = Double.parseDouble(tempStringSplit[5]);  
                        break;

                    case 7:

                        chr7[probeChr7++][col] = Double.parseDouble(tempStringSplit[5]);  
                        break;

                    case 8:

                        chr8[probeChr8++][col] = Double.parseDouble(tempStringSplit[5]);  
                        break;                   

                    case 9:

                        chr9[probeChr9++][col] = Double.parseDouble(tempStringSplit[5]);  
                        break;

                    case 10:

                        chr10[probeChr10++][col] = Double.parseDouble(tempStringSplit[5]);  
                        break;

                    case 11:

                        chr11[probeChr11++][col] = Double.parseDouble(tempStringSplit[5]);  
                        break;                        

                    case 12:

                        chr12[probeChr12++][col] = Double.parseDouble(tempStringSplit[5]);  
                        break;

                    case 13:

                        chr13[probeChr13++][col] = Double.parseDouble(tempStringSplit[5]);  
                        break;

                    case 14:

                        chr14[probeChr14++][col] = Double.parseDouble(tempStringSplit[5]);  
                        break;                        

                    case 15:

                        chr15[probeChr15++][col] = Double.parseDouble(tempStringSplit[5]);  
                        break;

                    case 16:

                        chr16[probeChr16++][col] = Double.parseDouble(tempStringSplit[5]);  
                        break;

                    case 17:

                        chr17[probeChr17++][col] = Double.parseDouble(tempStringSplit[5]); 
                        break;                        

                    case 18:

                        chr18[probeChr18++][col] = Double.parseDouble(tempStringSplit[5]); 
                        break;

                    case 19:

                        chr19[probeChr19++][col] = Double.parseDouble(tempStringSplit[5]); 
                        break;

                    case 20:

                        chr20[probeChr20++][col] = Double.parseDouble(tempStringSplit[5]); 
                        break;

                    case 21:

                        chr21[probeChr21++][col] = Double.parseDouble(tempStringSplit[5]); 
                        break;                        

                    case 22:

                        chr22[probeChr22++][col] = Double.parseDouble(tempStringSplit[5]); 
                        break;
                } 

            }
        col++;
        
        }
        
        BufferedWriter writeOut = null;
        
        int probeID = 0;
     
        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr1"));
        System.out.println("Writing segmented copy number matrix for chromosome1... ");
        Logger.logging("Writing segmented copy number matrix for chromosome1... ");

        int j;
        for (j=0; j<Chromosome.LENCHR1; j++)
        {                
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";                   
            probeID++;            
            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr1[j][k] + "\t"; 
            }      
            writeOut.write(entry);                
            writeOut.newLine();                
        }
        writeOut.close();

        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr2"));
        System.out.println("Writing segmented copy number matrix for chromosome2... ");
        Logger.logging("Writing segmented copy number matrix for chromosome2... ");

        for (j=0; j<Chromosome.LENCHR2; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";                 
            probeID++;            
            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr2[j][k] + "\t"; 
            }

            writeOut.write(entry);
            writeOut.newLine();                
        }
         writeOut.close();
        
        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr3"));
        System.out.println("Writing segmented copy number matrix for chromosome3... ");
        Logger.logging("Writing segmented copy number matrix for chromosome3... ");

        for (j=0; j<Chromosome.LENCHR3; j++)
        {
           String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";   

            probeID++;

            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr3[j][k] + "\t"; 
            }            
            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();
        
        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr4"));
        System.out.println("Writing segmented copy number matrix for chromosome4... ");
        Logger.logging("Writing segmented copy number matrix for chromosome4... ");

        for (j=0; j<Chromosome.LENCHR4; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";   

            probeID++;

            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr4[j][k] + "\t"; 
            }

            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();
        
        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr5"));
        System.out.println("Writing segmented copy number matrix for chromosome5... ");
        Logger.logging("Writing segmented copy number matrix for chromosome5... ");

        for (j=0; j<Chromosome.LENCHR5; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";   

            probeID++;

            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr5[j][k] + "\t"; 
            }

            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();
        
        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr6"));
        System.out.println("Writing segmented copy number matrix for chromosome6... ");
        Logger.logging("Writing segmented copy number matrix for chromosome6... ");

        for ( j=0; j<Chromosome.LENCHR6; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";   

            probeID++;

            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr6[j][k] + "\t"; 
            }

            writeOut.write(entry);
            writeOut.newLine();

        }
        writeOut.close();

        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr7"));
        System.out.println("Writing segmented copy number matrix for chromosome7... ");
        Logger.logging("Writing segmented copy number matrix for chromosome7... ");

        for (j=0; j<Chromosome.LENCHR7; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";   

            probeID++;

            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr7[j][k] + "\t"; 
            }

            writeOut.write(entry);
            writeOut.newLine();

        }
        writeOut.close();

        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr8"));
        System.out.println("Writing segmented copy number matrix for chromosome8... ");
        Logger.logging("Writing segmented copy number matrix for chromosome8... ");

        for (j=0; j<Chromosome.LENCHR8; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";   

            probeID++;

            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr8[j][k] + "\t"; 
            }

            writeOut.write(entry);
            writeOut.newLine();

        }
        writeOut.close();

        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr9"));
        System.out.println("Writing segmented copy number matrix for chromosome9... ");
        Logger.logging("Writing segmented copy number matrix for chromosome9... ");

        for (j=0; j<Chromosome.LENCHR9; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";                   
            probeID++;            
            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr9[j][k] + "\t"; 
            }            
            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();

        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr10"));
        System.out.println("Writing segmented copy number matrix for chromosome10... ");
        Logger.logging("Writing segmented copy number matrix for chromosome10... ");

        for (j=0; j<Chromosome.LENCHR10; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";                 
            probeID++;            
            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr10[j][k] + "\t"; 
            }            
            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();
       
        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr11"));
        System.out.println("Writing segmented copy number matrix for chromosome11... ");
        Logger.logging("Writing segmented copy number matrix for chromosome11... ");

        for (j=0; j<Chromosome.LENCHR11; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";                 
            probeID++;            
            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr11[j][k] + "\t"; 
            }            
            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();
        
        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr12"));
        System.out.println("Writing segmented copy number matrix for chromosome12... ");
        Logger.logging("Writing segmented copy number matrix for chromosome12... ");

        for (j=0; j<Chromosome.LENCHR12; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";                   
            probeID++;            
            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr12[j][k] + "\t"; 
            }            
            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();

        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr13"));
        System.out.println("Writing segmented copy number matrix for chromosome13... ");
        Logger.logging("Writing segmented copy number matrix for chromosome13... ");

        for (j=0; j<Chromosome.LENCHR13; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";                  
            probeID++;            
            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr13[j][k] + "\t"; 
            }            
            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();
        
        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr14"));
        System.out.println("Writing segmented copy number matrix for chromosome14... ");
        Logger.logging("Writing segmented copy number matrix for chromosome14... ");

        for (j=0; j<Chromosome.LENCHR14; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";                  
            probeID++;            
            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr14[j][k] + "\t"; 
            }            
            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();

        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr15"));
        System.out.println("Writing segmented copy number matrix for chromosome15... ");
        Logger.logging("Writing segmented copy number matrix for chromosome15... ");

        for (j=0; j<Chromosome.LENCHR15; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";                  
            probeID++;            
            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr15[j][k] + "\t"; 
            }            
            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();
        
        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr16")); 
        System.out.println("Writing segmented copy number matrix for chromosome16... ");
        Logger.logging("Writing segmented copy number matrix for chromosome16... ");
        for (j=0; j<Chromosome.LENCHR16; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";                  
            probeID++;            
            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr16[j][k] + "\t"; 
            }            
            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();
        
        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr17"));
        System.out.println("Writing segmented copy number matrix for chromosome17... ");
        Logger.logging("Writing segmented copy number matrix for chromosome17... ");

        for (j=0; j<Chromosome.LENCHR17; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";                  
            probeID++;            
            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr17[j][k] + "\t"; 
            }            
            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();
               
        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr18"));
        System.out.println("Writing segmented copy number matrix for chromosome18... ");
        Logger.logging("Writing segmented copy number matrix for chromosome18... ");

        for (j=0; j<Chromosome.LENCHR18; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";                  
            probeID++;            
            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr18[j][k] + "\t"; 
            }            
            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();
        
        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr19"));
        System.out.println("Writing segmented copy number matrix for chromosome19... ");
        Logger.logging("Writing segmented copy number matrix for chromosome19... ");

        for (j=0; j<Chromosome.LENCHR19; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";                 
            probeID++;            
            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr19[j][k] + "\t"; 
            }
            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();
        
        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr20"));
        System.out.println("Writing segmented copy number matrix for chromosome20... ");
        Logger.logging("Writing segmented copy number matrix for chromosome20... ");

        for (j=0; j<Chromosome.LENCHR20; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";                
            probeID++;            
            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr20[j][k] + "\t"; 
            }            
            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();
      
        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr21"));
        System.out.println("Writing segmented copy number matrix for chromosome21... ");
        Logger.logging("Writing segmented copy number matrix for chromosome21... ");

        for (j=0; j<Chromosome.LENCHR21; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";                  
            probeID++;            
            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr21[j][k] + "\t"; 
            }            
            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();

        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "Chr22"));
        System.out.println("Writing segmented copy number matrix for chromosome22... ");
        Logger.logging("Writing segmented copy number matrix for chromosome22... ");

        for (j=0; j<Chromosome.LENCHR22; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";                  
            probeID++;            
            for (int k=0; k<numSample; k++)
            {
                entry = entry + chr22[j][k] + "\t"; 
            }            
            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();

        probeID = 0;
        
        writeOut = new BufferedWriter(new FileWriter(saicInputDir + "GenomeWide"));
        System.out.println("Writing segmented copy number matrix for the whole genome... ");
        Logger.logging("Writing segmented copy number matrix for the whole genome... ");

        for (j=0; j<Chromosome.NUMPROBE; j++)
        {
            String entry = chrID[probeID] + "\t" + probeName[probeID] + "\t" + probeLoci[probeID]+ "\t";                  
            probeID++;            
            for (int k=0; k<numSample; k++)
            {
                entry = entry + wholeGenome[j][k] + "\t"; 
            }            
            writeOut.write(entry);
            writeOut.newLine();                
        }
        writeOut.close();

        chr1 = null;
        chr2 = null;
        chr3 = null;
        chr4 = null;
        chr5 = null;
        chr6 = null;
        chr7 = null;
        chr8 = null;
        chr9 = null;
        chr10 = null;
        chr11 = null;
        chr12 = null;
        chr13 = null;
        chr14 = null;
        chr15 = null;
        chr16 = null;
        chr17 = null;
        chr18 = null;
        chr19 = null;
        chr20 = null;
        chr21 = null;
        chr22 = null; 
        return;
    }

}
    
    
