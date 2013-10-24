/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.vt.cbil.ui;

import edu.vt.cbil.bacom.CorrectNormalContam;
import edu.vt.cbil.exceptions.EmptyStringException;
import edu.vt.cbil.exceptions.FileDoesNotExistException;
import edu.vt.cbil.exceptions.InvalidFormatException;
import edu.vt.cbil.exceptions.InvalidInputFileException;
import edu.vt.cbil.exceptions.SampleSizeDoesNotMatchException;
import edu.vt.cbil.saic.AisaicPlot;
import static edu.vt.cbil.saic.AisaicPlot.createAisaicPlotObjectFromDataFiles;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.Parser;
import org.apache.commons.cli.PosixParser;
import edu.vt.cbil.saic.DetectSCAs;
import edu.vt.cbil.util.BacomParameters;
import edu.vt.cbil.util.DataConversion;
import edu.vt.cbil.util.GeneSet;
import edu.vt.cbil.util.Logger;
import edu.vt.cbil.util.SaicParameters;

/**
 *
 * @author bella
 */
public class AisaicCmdLineUI {  
    
    private static  Options createOptions(){
        Options options = new Options();
        options.addOption("t", "analysis type", true, "choose the analysis type (required). Values: 0 -- BACOM; 1 -- SAIC; 2 -- AISAIC");
        options.addOption("mb", "BACOM analysis mode", true, "analysis mode for BACOM (required when BACOM or AISAIC is chosen). Values: 0 -- single sample analysis; 1 -- multiple samples analysis");
        options.addOption("ns", "normal samples", true, "the path of normal samples (required for BACOM or AISAIC analysis). For single sample analysis, it is the directory of the normal sample; "
                + "for multiple-sample analysis, it is the directory containing all the normal samples AND a text file specifying all the normal sample file names to be analyzed.");
        options.addOption("ts", "tumor samples", true, "the directory of tumor samples (required for BACOM or AISAIC analysis). For single sample analysis is chosen, it is the directory of the tumor sample;"
                + "for multiple-sample analysis, it is the directory containing all the tumor samples AND a text file specifying all the tumor sample file names to be analyzed.");
        options.addOption("cdf", "cdf file", true, "the directory of the cdf file (required for BACOM or AISAIC analysis).");
        options.addOption("af", "annotation file", true, "the directory of the annotation file (required for BACOM or AISAIC analysis).");
        options.addOption("rb", "result folder for BACOM/AISAIC", true, "the path of the result folder for BACOM/AISAIC (required for BACOM or AISAIC). If the specified path doesn't have folder named 'result', the program will create one.");
        options.addOption("gi", "genomic region of interest", true, "the genomic regions of interest (optional, provided only when user wants to know the deletion type of this region). The value format: 'chr:startLoci-endLoci'");
        options.addOption("ms", "SAIC analysis mode", true, "analysis mode for SAIC (required when SAIC or AISAIC is chosen). Values: 0 -- single chromosome analysis; 1 -- multiple chromosomes analysis; 2 -- genome-wide analysis");
        options.addOption("sd", "segmented CNA data", true, "the directory of segmented copy number data matrix (required only for SAIC). The data format is a two dimensional matrix of which each column is the segmented copy number profile of a sample.");
        options.addOption("rs", "result folder for SAIC", true, "the directory of the result folder for SAIC (required only when SAIC is chosen). If the specified path doesn't have folder named 'result', the program will create one.");
        options.addOption("gr", "gene reference file", true, "the directory of the reference gene. This is used to query the genes covered by detected SCAs.");
        options.addOption("ta", "amp threshold", true, "the threshold for detecting CNA amplification probes (optional). The default value is 0.322.");
        options.addOption("td", "del threshold", true, "the threshold for detecting CNA deletion probes (optional). The default value is -0.415.");
        options.addOption("np", "number of permutations", true, "the number of permutations (optional). The default value is 1000.");
        options.addOption("tr", "correlation coefficient threshold", true, "the threshold of correlation coefficient between consecutive probes, used to construct CNA units (optional). The default value is 0.9.");
        options.addOption("ip", "iterative permutation", true, "a binary value. Values: ip=0 -- without iterative permutation; ip=1 -- with iterative permutation. Default value is ip=1");
        options.addOption("qs", "quick SAIC", true, "a binary value. Values: q=1 -- with downsampling to speed up the analysis; q=0 -- without downsampling. Default value is q=1");
        options.addOption("h", "help", true, "print help information");
        return options;
        
    }
   
    private static void showHelp(Options options){
        HelpFormatter h = new HelpFormatter();
        h.printHelp("help", options);
        System.exit(-1);
    }
    
    public static void main(String[] args){

        Options options = createOptions();
    //    showHelp(options);
        try{
            CommandLineParser  cmdParser = new GnuParser();
            CommandLine cmd = cmdParser.parse(options, args);
            
            if (cmd.hasOption("h")){
                showHelp(options);
                System.exit(-1);
            }
            
            if (!cmd.hasOption("t")){
                System.out.println("You need to specify the analysis you want");
                showHelp(options);
            }
            else{
                if (cmd.getOptionValue("t").equals("0")){
                    System.out.println("BACOM analysis is chosen");
                    startBACOM(cmd, options);
                }
                    
                else if (cmd.getOptionValue("t").equals("1")){
                    System.out.println("SAIC analysis is chosen");
                    startSAIC(cmd, options);  
                }
                
                else if (cmd.getOptionValue("t").equals("2")){
                    System.out.println("AISAIC analysis is chosen");
                    startAISAIC(cmd, options);
                }
                
                else {
                    System.out.println("Wrong analysis type");
                    showHelp(options);
                }      
            }      
        }catch (Exception e){
            e.printStackTrace();
            showHelp(options);             
        }
 
    }
    
    public static void startBACOM(CommandLine cmd, Options options){
       
        BacomParameters bacomPara = new BacomParameters();
        if (!cmd.hasOption("mb")){
            System.out.println("Please specify the BACOM analysis mode");
        }
        else{
            if (cmd.getOptionValue("mb").equals("0")){
                bacomPara.setAnalysisMode("0");
                System.out.println("Single sample BACOM analysis is chosen" );
            }
            else if (cmd.getOptionValue("mb").equals("1")){
                System.out.println("Multiple samples BACOM analysis is chosen");
                bacomPara.setAnalysisMode("1");      
            }
            else{
                System.out.println("Please specify the right BACOM analysis mode, could only be '0' or '1'");
                System.exit(-1);
            }
            if (!cmd.hasOption("ns")){
                System.out.println("Please specify a normal file");
                System.exit(-1);
            }
            else{                                  
                bacomPara.setInputNormalFile(cmd.getOptionValue("ns"));  
                try{
                    bacomPara.setNormalSamples();
                }
                catch (InvalidInputFileException e){
                    System.out.println("The input normal sample should be a .CEL file for single sample analysis");  
                }
                catch (IOException e1){
                    System.out.println("File not found/Error reading the file!" );              
                }
                catch (FileDoesNotExistException e2){
                    System.out.println("File does not exist, please check again." );
                }
                System.out.println("The selected normal file is: " + bacomPara.getInputNormalFile());
            }
            if (!cmd.hasOption("ts")){ 
                System.out.println("Please specify tumor file");
                System.exit(-1);
            }
            else{                                
                bacomPara.setInputTumorFile(cmd.getOptionValue("ts")); 
                try{
                    bacomPara.setTumorSamples();
                }catch (SampleSizeDoesNotMatchException e){
                    System.out.println("Sample size does not MATCH! Plese check your input file" );
                    System.exit(-1);
                }
                catch (InvalidInputFileException e2){
                    System.out.println("The input tumor sample should be a .CEL file for single sample analysis" );    
                    System.exit(-1);
                }
                catch (IOException e3){
                    System.out.println("File not found/Error reading the file!");  
                    System.exit(-1);
                }
                catch (FileDoesNotExistException e4){
                    System.out.println("File does not exist in the specified folder, please check again.");
                    System.exit(-1);
                }
                System.out.println("The selected tumor file is: " + bacomPara.getInputNormalFile());
            }
            if (!cmd.hasOption("cdf")){
                System.out.println("Please specify the CDF file" );
                System.exit(-1);
            }
            else{  
                try{
                    bacomPara.setInputCdfFile(cmd.getOptionValue("cdf"));  
                }catch (FileDoesNotExistException e){
                    System.out.println(e);
                    System.exit(-1);
                }
                System.out.println("The selected CDF file is: " + bacomPara.getInputCdfFile());
            }
            if (!cmd.hasOption("af")){
                System.out.println("Please specify the annotation file" );
                System.exit(-1);
            }
            else{ 
                try{
                    bacomPara.setAnnotationFile(cmd.getOptionValue("af")); 
                }catch (FileDoesNotExistException e){
                    System.out.println(e);
                    System.exit(-1);
                } 
                System.out.println("The selected annotation file is: " + bacomPara.getInputCdfFile());
            }

            if (!cmd.hasOption("rb")){
                System.out.println("Please specify the directory where to store the BACOM results");
                System.exit(-1);
            }
            else{
                File resultDir = new File(cmd.getOptionValue("rb") + File.separator + "results");
                if (!resultDir.exists())
                    resultDir.mkdir();
                bacomPara.setResultFolder(cmd.getOptionValue("rb")+ File.separator + "results" + File.separator);
                System.out.println("The analysis result will be stored in 'result' folder under " + bacomPara.getResultFolder());                
                File loggerFile = new File(bacomPara.getResultFolder() + File.separator + "logFile.txt");
                if (loggerFile != null && loggerFile.exists())
                    loggerFile.delete();   
                try{
                    Logger.setLogger(loggerFile.getPath());
                }catch (IOException e){
                    System.out.println(e);
                }
            }
            if (cmd.hasOption("gi")){
                System.out.println("The interested genomic region is: " + cmd.getOptionValue("gi") );
                try{
                    System.out.println("The interested genomic region is " + cmd.getOptionValue("gi"));
                    bacomPara.setGeneInterested(cmd.getOptionValue("gi"));
                }catch (InvalidFormatException e){
                    System.out.println("Invalid format for interested genomic region, the right format should be like: '10:123456-456789'");
                    System.exit(-1);
                }
                catch(EmptyStringException e2){
                    System.out.println("User choose -gi option, but didn't specify the interested genomic region,please check.");
                    System.exit(-1);
                }
            }
            else{
                System.out.println("User didn't specify the interested genomic region" );
                bacomPara.setGeneInterestedNull();
            }

            if (bacomPara.getAnalysisMode().equals("0")){
                Long startTime = System.currentTimeMillis();
                System.out.println("Start estimating the normal tissue contamination for tumor sample: " + bacomPara.getSampleIds().get(0));
                CorrectNormalContam singleCorrectNormalContam = new CorrectNormalContam(bacomPara, 0, 0);                
                ForkJoinPool bacomForkJoinPool = new ForkJoinPool();    
                bacomForkJoinPool.invoke(singleCorrectNormalContam);
                System.out.println("The time needed to analysis this sample is: " + (System.currentTimeMillis() - startTime)/1000 + " seconds");
            }
            else if (bacomPara.getAnalysisMode().equals("1")){
                
                int sampleSize = bacomPara.getSampleSize();
                int numSampleConcur = Runtime.getRuntime().availableProcessors();
                if (numSampleConcur > 4)
                    numSampleConcur = numSampleConcur -2;
                else if (numSampleConcur <= 2)
                    numSampleConcur = 1;
                if (sampleSize < numSampleConcur)
                {
                    numSampleConcur = sampleSize;
                }
                for (int j=0; j<sampleSize/numSampleConcur; j++)
                { 
                    System.out.println("Start parallelly estimating the normal tissue contamination for tumor samples: " );   
                    for (int k=numSampleConcur*j; k<numSampleConcur*(j+1); k++){
                        System.out.println(bacomPara.getSampleIds().get(k));
                    }
                    Long startTime = System.currentTimeMillis();
                    ForkJoinPool bacomForkJoinPoolMultiple = new ForkJoinPool();                        
                    CorrectNormalContam multiCorrectNormalContam = new CorrectNormalContam(bacomPara,numSampleConcur*j, numSampleConcur*(j+1)-1);               
                    bacomForkJoinPoolMultiple.invoke(multiCorrectNormalContam);
                    System.out.println("The time needed for analyzing these " + numSampleConcur + " samples is: " + (System.currentTimeMillis()-startTime)/1000 + " seconds.");
                }                            
                int remain = sampleSize%numSampleConcur; 
                System.out.println("Start estimating the normal tissue contamination for remaining tumor samples: " );  
                for (int k=sampleSize-remain; k<sampleSize-1; k++){
                        System.out.println(bacomPara.getSampleIds().get(k));
                    }
                ForkJoinPool bacomForkJoinPoolMultiple = new ForkJoinPool();                        
                CorrectNormalContam multiCorrectNormalContam = new CorrectNormalContam(bacomPara, sampleSize-remain, sampleSize-1);                    
                bacomForkJoinPoolMultiple.invoke(multiCorrectNormalContam);             
            }   
        } 
        System.out.println("Done with the BACOM analysis");
                
    }
    
   public static void startSAIC(CommandLine cmd, Options options){
       SaicParameters saicPara = new SaicParameters();
       if (!cmd.hasOption("ms")){
           System.out.println("Please specify the analysis mode for SAIC. Single chromosome analysis = 0; Multiple chromosomes analysis=1;Genome-wide analysis=2");          
           System.exit(-1);
       }
       else{
           if (cmd.getOptionValue("ms").equals("0")){
               System.out.println("User choose single chromosome analysis");
               saicPara.setAnalysisMode(0);
           }
           else if (cmd.getOptionValue("ms").equals("1")){
               System.out.println("User choose multiple chromosomes analysis");
               saicPara.setAnalysisMode(1);
           }
           else if (cmd.getOptionValue("ms").equals("2")){
               System.out.println("User choose genome-wide analysis");
               saicPara.setAnalysisMode(2);
           }
           else{
               System.out.println("Wrong analysis mode input, please try again. Valid value should be: 0, 1, 2");
               System.exit(-1);
           }
           if (!cmd.hasOption("sd")){
               System.out.println("Please specify the input segmented copy number data matrix");
               System.exit(-1);
           }
           else{
               if(saicPara.getAnalysisMode() == 1)
                   saicPara.setInfileDir(cmd.getOptionValue("sd"));
               else if (saicPara.getAnalysisMode() == 0 || saicPara.getAnalysisMode() == 2)
                   saicPara.setInfile(cmd.getOptionValue("sd"));
           }
           if (!cmd.hasOption("rs")){
               System.out.println("Please specify the path where the SCAs detection results would be stored");
               System.exit(-1);
           }
           else{
               File resultDir = new File(cmd.getOptionValue("rs") + File.separator + "SAICResults");
               if (!resultDir.exists())
                   resultDir.mkdir();
               saicPara.setOutfileDir(resultDir.getAbsolutePath());
               File loggerFile = new File(cmd.getOptionValue("rs")+ File.separator + "logFile.txt");
                if (loggerFile != null && loggerFile.exists())
                    loggerFile.delete();   
                try{
                    Logger.setLogger(loggerFile.getAbsolutePath());
                }catch (IOException e){
                    System.out.println(e);
                }                
           }
           if (!cmd.hasOption("gr")){
               System.out.println("Please choose the gene reference file");
               System.exit(-1);
           }
           else{
               try{
                    saicPara.setGeneRefDir(cmd.getOptionValue("gr"));
               }catch(FileDoesNotExistException e){
                   System.out.println("Can not find the specified gene reference file, please check again.");
                   System.exit(-1);
               }
               GeneSet.readGeneInfo(saicPara.getGeneRefDir());
               
           }                      
           if (!cmd.hasOption("ta")){
               System.out.println("User didn't specify the threshold for copy number amplification detection, the software will use the default value");
           }
           else{ 
               double ta = 0.0;
               try {                    
                    ta = Double.parseDouble(cmd.getOptionValue("ta"));
                    if (ta<0)
                    {
                        System.out.println("Invalid value for copy number deletion threshold. Please choose a positive float number");                     
                        System.exit(-1);
                    }
                }catch (NumberFormatException e)
                {
                    System.out.println("Invalid value for copy number deletion threshold. Please choose a positive float number");
                    System.exit(-1);
                }
               saicPara.setTa(ta);
               System.out.println("The threshold for copy number amplification detection is: " + cmd.getOptionValue("ta"));
           }
           if (!cmd.hasOption("td")){
               System.out.println("User didn't specify the threshold for copy number deletion detection, the software will use the default value");
           }
           else{ 
               double td = -0.0;
               try {
                    td = Double.parseDouble(cmd.getOptionValue("td"));
                    td = 0- td;
                    if (td>0)
                    {
                        System.out.println("Invalid value for copy number deletion threshold. Please choose a negative float number");                     
                        System.exit(-1);
                    }
                }catch (NumberFormatException e)
                {
                    System.out.println("Invalid value for copy number deletion threshold. Please choose a negative float number" );
                    System.exit(-1);
                }
               saicPara.setTd(td);
               System.out.println("The threshold for copy number deletion detection is: " + saicPara.getTd());
           }
           
           if (!cmd.hasOption("np")){
               System.out.println("User didn't specify the number of permutation, the software will use the default value");
           }
           else{ 
               int np = 0;
               try {
                    np = Integer.parseInt(cmd.getOptionValue("np"));
                    if (np<0)
                    {
                        System.out.println("Invalid value for number of permutation. Please choose in a positive integer.");                     
                        System.exit(-1);
                    }
                }catch (NumberFormatException e)
                {
                    System.out.println("Invalid value for number of permutation. Please choose in a positive integer." );
                    System.exit(-1);
                }
               saicPara.setNp(np);
               System.out.println("The number of permutation is: " + saicPara.getNp());
           }
           
           if (!cmd.hasOption("tr")){
               System.out.println("User didn't specify the threshold for correltion coefficient between consecutive probes, the software will use the default value");
           }
           else{ 
               double tr = 0.0;
               try {
                    tr = Double.parseDouble(cmd.getOptionValue("tr"));
                    if (tr<0 || tr>1)
                    {
                        System.out.println("Invalid value for correlation coefficient threshold. Please type in a float number between 0 and 1 " );                     
                        System.exit(-1);
                    }
                }catch (NumberFormatException e)
                {
                    System.out.println("Invalid value for correlation coefficient threshold. Please type in a float number between 0 and 1 ");                     
                    System.exit(-1);
                }
               saicPara.setTci(tr);
               System.out.println("The threshold for correlation coefficient is: " + saicPara.getTci());
           }
           if (!cmd.hasOption("ip")){
               System.out.println("User didn't specify whether to use SCA-exclusive permutation, the program will use this scheme as default" );
           }
           else{
               if (cmd.getOptionValue("ip").equals("0")){
                   saicPara.setIterFlag(false);
                   System.out.println("User doesn't select use SCA-exclusive permutation scheme.");
               }
               else{
                   saicPara.setIterFlag(true);
                   System.out.println("User selects SCA-exclusive permutation scheme.");
               }
           }
           if (!cmd.hasOption("qs")){
               System.out.println("User didn't specify whether to use quick SCAs detection, the program will adopt this as default" );
           }
           else{
               if (cmd.getOptionValue("qs").equals("0")){
                   saicPara.setQuickSAIC(false);
                   System.out.println("User doesn't select quick SAIC analysis scheme.");
               }
               else{
                   saicPara.setQuickSAIC(true);
                   System.out.println("User selects quick SAIC analysis scheme.");
               }
           }
           if (saicPara.getAnalysisMode() == 0 || saicPara.getAnalysisMode() == 2){
                String infiles[] = new String[1];
                String outfiles[] = new String[1];
                infiles[0] = saicPara.getInfile();
                if (saicPara.getAnalysisMode() == 0)
                    outfiles[0] = saicPara.getOutfileDir() + File.separator + "Result";
                if (saicPara.getAnalysisMode() == 2)
                    outfiles[0] = saicPara.getOutfileDir() + File.separator + "GenomeWideResult";
                Long startTime = System.currentTimeMillis();
                DetectSCAs detectSCAs = new DetectSCAs(saicPara, infiles, outfiles);                
                ForkJoinPool detectSCAForkJoinPool = new ForkJoinPool();                                       
                detectSCAForkJoinPool.invoke(detectSCAs); 
                if (saicPara.getAnalysisMode() == 0)
                    System.out.println("The analyzing time for this chromosome is " + (System.currentTimeMillis()-startTime)/1000 + " seconds");
                if (saicPara.getAnalysisMode() == 2)
                    System.out.println("The analyzing time for the whole genome is " + (System.currentTimeMillis()-startTime)/1000 + " seconds");
                detectSCAs = null;
                detectSCAForkJoinPool = null;
           }
           else if (saicPara.getAnalysisMode() == 1){
                String infiles[] = new String[1];
                String outfiles[] = new String[1];
                Long startTime = 0L;
               for (int chr =22; chr>0; chr--) 
                {
                    infiles[0] = saicPara.getInfileDir()+ File.separator + "Chr" + chr;
                    outfiles[0] = saicPara.getOutfileDir()+ File.separator + "ResultChr" + chr;                                                         
                    DetectSCAs detectSCAs = new DetectSCAs(saicPara, infiles, outfiles);                                   
                    ForkJoinPool detectSCAForkJoinPool = new ForkJoinPool();                    
                    startTime = System.currentTimeMillis();                    
                    detectSCAForkJoinPool.invoke(detectSCAs);                                        
                    System.out.println("The analyzing time for chromosome" + chr + " is " + (System.currentTimeMillis()-startTime)/1000 + " seconds");
                    detectSCAs = null;
                    detectSCAForkJoinPool = null;
                }       
           }
           // Plot SAIC result figures

            double[] cn;
            String[] chrID;
            int[] loc;
            String fileName;
            double[] meanVal;
            boolean[] isSignificant;

            final String DATA_PATH = saicPara.getOutfileDir() + File.separator;
            final String AMP_FILE_SUFFIX = "_ampAvgCNA";
            final String DEL_FILE_SUFFIX = "_delAvgCNA";

            final String FIG_PATH = DATA_PATH + "figures" + File.separator;

            File figureDir = new File(FIG_PATH);

            // if the directory does not exist, create it
            if (!figureDir.exists()) {
                figureDir.mkdir();  
            }
            if (saicPara.getAnalysisMode() ==0 ){
                final String FILE_PREFIX = "Result";
                String ampFile = DATA_PATH + FILE_PREFIX +  AMP_FILE_SUFFIX;
                String delFile = DATA_PATH + FILE_PREFIX + DEL_FILE_SUFFIX;
                String figureFileName = FIG_PATH + "AISAIC_singleChr.png";
                String chr0 = "";
                AisaicPlot p = createAisaicPlotObjectFromDataFiles(chr0,
                        ampFile, delFile);
                p.plot();
                p.saveChartAsPNG(figureFileName);                
            }
            else if ( saicPara.getAnalysisMode() == 2){
                final String FILE_PREFIX_GENOME = "ResultGenomeWide";
                String ampFile = DATA_PATH + FILE_PREFIX_GENOME +  AMP_FILE_SUFFIX;
                String delFile = DATA_PATH + FILE_PREFIX_GENOME + DEL_FILE_SUFFIX;
                String figureFileName = FIG_PATH + "AISAIC_genome.png";
                String genome = "genome";
                AisaicPlot p = createAisaicPlotObjectFromDataFiles(genome,
                        ampFile, delFile);
                p.plot();
                p.saveChartAsPNG(figureFileName);
            }

            else if(saicPara.getAnalysisMode() == 1){
                String[] chrArray = {"1", "2", "3", "4", "5", "6", "7", "8", 
                "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", 
                "19", "20", "21", "22"};
                for (String chr: chrArray) {
                    final String FILE_PREFIX_CHR = "ResultChr";
                    String ampFile = DATA_PATH + FILE_PREFIX_CHR + chr + AMP_FILE_SUFFIX;
                    String delFile = DATA_PATH + FILE_PREFIX_CHR + chr + DEL_FILE_SUFFIX;
                    String figureFileName = FIG_PATH + "AISAIC_chr_" + chr + ".png";
                    AisaicPlot p = createAisaicPlotObjectFromDataFiles(chr,
                            ampFile, delFile);
                    p.plot();
                    p.saveChartAsPNG(figureFileName);
                }
            }

       }
       
    }
        
   public static void startAISAIC(CommandLine cmd, Options options){
       
       BacomParameters bacomPara = new BacomParameters();
       SaicParameters saicPara = new SaicParameters();
       File loggerFile;
       if (!cmd.hasOption("mb")){
            System.out.println("User didn't specify the BACOM analysis mode, the program will use 'multiple samples analysis by default'");
            bacomPara.setAnalysisMode("1");
       }
        else{
            if (cmd.getOptionValue("mb").equals("0")){
                System.out.println("Please choose multiple samples analysis for your convenience" );
                System.exit(-1);
            }
            else if (cmd.getOptionValue("mb").equals("1")){
                System.out.println("Multiple samples BACOM analysis is chosen");
                bacomPara.setAnalysisMode("1");      
            }
            else{
               System.out.println("Please set valid BACOM analysis mode.");
               System.exit(-1);
            }
       }
            if (!cmd.hasOption("ns")){
                System.out.println("Please specify a normal file");
                System.exit(-1);
            }
            else{                                  
                bacomPara.setInputNormalFile(cmd.getOptionValue("ns"));  
                try{
                    bacomPara.setNormalSamples();
                }
                catch (InvalidInputFileException e){
                    System.out.println("The input normal sample should be a .CEL file for single sample analysis"); 
                    System.exit(-1);
                }
                catch (IOException e1){
                    System.out.println("File not found/Error reading the file!");     
                    System.exit(-1);
                }
                catch (FileDoesNotExistException e2){
                    System.out.println("File does not exist in the specified folder, please check again");
                    System.exit(-1);
                }
                System.out.println("The selected normal file is: " + bacomPara.getInputNormalFile());
            }
            if (!cmd.hasOption("ts")){ 
                System.out.println("Please specify tumor file" );
                System.exit(-1);
            }
            else{                                
                bacomPara.setInputTumorFile(cmd.getOptionValue("ts")); 
                try{
                    bacomPara.setTumorSamples();
                }catch (SampleSizeDoesNotMatchException e){
                    System.out.println("Sample size does not MATCH! Plese check your input file" );
                    System.exit(-1);
                }
                catch (InvalidInputFileException e2){
                    System.out.println("The input tumor sample should be a .CEL file for single sample analysis" );  
                    System.exit(-1);
                }
                catch (IOException e3){
                    System.out.println("File not found/Error reading the file!" );  
                    System.exit(-1);
                }
                catch (FileDoesNotExistException e2){
                    System.out.println("File does not exist in the specified folder, please check again");
                    System.exit(-1);
                }
                System.out.println("The selected tumor file is: " + bacomPara.getInputNormalFile());
            }
            if (!cmd.hasOption("cdf")){
                System.out.println("Please specify the CDF file");
                System.exit(-1);
            }
            else{  
                try{
                    bacomPara.setInputCdfFile(cmd.getOptionValue("cdf"));  
                }catch (FileDoesNotExistException e){
                    System.out.println(e);
                }
                System.out.println("The selected CDF file is: " + bacomPara.getInputCdfFile());
            }
            if (!cmd.hasOption("af")){
                System.out.println("Please specify the annotation file" );
                System.exit(-1);
            }
            else{ 
                try{
                    bacomPara.setAnnotationFile(cmd.getOptionValue("af")); 
                }catch (FileDoesNotExistException e){
                    System.out.println(e);
                } 
                System.out.println("The selected annotation file is: " + bacomPara.getInputCdfFile());
            }

            if (!cmd.hasOption("rb")){
                System.out.println("Please specify the directory where to store the BACOM results");
                System.exit(-1);
            }
            else{                
                File resultDir = new File(cmd.getOptionValue("rb") + File.separator + "results");
                if(!resultDir.exists()){
                    resultDir.mkdir();
                } 
                bacomPara.setResultFolder(resultDir.getAbsolutePath());
                loggerFile = new File(bacomPara.getResultFolder() + File.separator + "logFile.txt");
                if (loggerFile != null && loggerFile.exists())
                    loggerFile.delete();   
                try{
                    Logger.setLogger(loggerFile.getPath());
                }catch (IOException e){
                    System.out.println(e);
                }
            }
            if (cmd.hasOption("gi")){
                System.out.println("You should be able to find the deletion type of this genomic region after analysis" );
                try{
                    bacomPara.setGeneInterested(cmd.getOptionValue("gi"));
                }catch (InvalidFormatException e){
                    System.out.println("Invalid format for interested genomic region, the right format should be like: '10:123456-456789'");
                    System.exit(-1);
                }
                catch(EmptyStringException e2){
                    System.out.println("User choose -gi option, but didn't specify the interested genomic region,please check.");
                    System.exit(-1);
                }
            }
            else{
                System.out.println("User doesn't have any interested genomic region" );
                bacomPara.setGeneInterestedNull();
            }

            if (bacomPara.getAnalysisMode().equals("1")){                
                int sampleSize = bacomPara.getSampleSize();
                int numSampleConcur = Runtime.getRuntime().availableProcessors();
                if (numSampleConcur > 4)
                    numSampleConcur = numSampleConcur -2;
                else if (numSampleConcur <= 2)
                    numSampleConcur = 1;
                if (sampleSize < numSampleConcur)
                {
                    numSampleConcur = sampleSize;
                }
                System.out.println("There are " + bacomPara.getSampleSize() + " samples in total");
                for (int j=0; j<sampleSize/numSampleConcur; j++)
                {
                    System.out.println("Start parallelly estimating the normal tissue contamination for tumor samples: " );   
                    for (int k=numSampleConcur*j; k<numSampleConcur*(j+1); k++){
                        System.out.println(bacomPara.getSampleIds().get(k));
                    }
                    System.out.println("This may take couple of minutes...");
                    Long startTime = System.currentTimeMillis();
                    ForkJoinPool bacomForkJoinPoolMultiple = new ForkJoinPool();                        
                    CorrectNormalContam multiCorrectNormalContam = new CorrectNormalContam(bacomPara,numSampleConcur*j, numSampleConcur*(j+1)-1);               
                    bacomForkJoinPoolMultiple.invoke(multiCorrectNormalContam);
                    System.out.println("The time needed for analyzing these " + numSampleConcur + " samples is: " + (System.currentTimeMillis()-startTime)/1000 + " seconds");
                }                            
                int remain = sampleSize%numSampleConcur; 
                System.out.println("Analyzing the remaining " + remain + " samples");
                ForkJoinPool bacomForkJoinPoolMultiple = new ForkJoinPool();                        
                CorrectNormalContam multiCorrectNormalContam = new CorrectNormalContam(bacomPara, sampleSize-remain, sampleSize-1);                    
                bacomForkJoinPoolMultiple.invoke(multiCorrectNormalContam);             
            }   
//        }
       
        //       Data format conversion from BACOM results to SAIC input               //
             
        saicPara.setBACOMResultDir(bacomPara.getResultFolder() + File.separator + "BACOMResults" + File.separator);
        
        File resultFolder = new File (bacomPara.getResultFolder() + File.separator + "SAICInput");
        if(!resultFolder.exists()){
            resultFolder.mkdir();
        }
        saicPara.setInfileDir(resultFolder.getAbsolutePath()+ File.separator);
        
        resultFolder = new File (bacomPara.getResultFolder() + File.separator + "SAICResults");
        if(!resultFolder.exists()){
            resultFolder.mkdir();
        }
        saicPara.setOutfileDir(resultFolder.getAbsolutePath() + File.separator);
        
        DataConversion convertForSAIC = new DataConversion(bacomPara.getAlpha(), bacomPara.getSampleIds(), saicPara.getBACOMResultDir(), saicPara.getInfileDir());             

         try{
            convertForSAIC.transform();
         }catch(IOException e)
         {
             System.out.println("Can not convert BACOM results to SAIC input because no BACOM results available.");
             System.exit(-1);            
         }  
         
         //       Start SAIC analysis       //
         
       if (!cmd.hasOption("ms")){
           System.out.println("Please specify the analysis mode for SAIC. Valid values are: Multiple chromosomes analysis=1;Genome-wide analysis=2" );          
           System.exit(-1);
       }
       else{
           if (cmd.getOptionValue("ms").equals("0")){
               System.out.println("Single chromosome analysis is not supported in AISAIC analysis" );
           }
           else if (cmd.getOptionValue("ms").equals("1")){
               System.out.println("User choose multiple chromosomes analysis");
               saicPara.setAnalysisMode(1);
           }
           else if (cmd.getOptionValue("ms").equals("2")){
               System.out.println("User choose genome-wide analysis" );
               saicPara.setAnalysisMode(2);
           }
           else{
               System.out.println("Please specify the analysis mode for SAIC analysis.");
           }
           if (cmd.hasOption("sd")){
               System.out.println("In AISAIC analysis, no segmented copy number data matrix needed.");
               System.exit(-1);
           }
           if (cmd.hasOption("rs")){
               System.out.println("In AISAIC analysis, no need to specify the result folder, the software will put the results in the same folder as BACOM's result." );
               //System.exit(-1);
           }
           if (!cmd.hasOption("gr")){
               System.out.println("Please choose the gene reference file" );
               System.exit(-1);
           }
           else{
               try{
                    saicPara.setGeneRefDir(cmd.getOptionValue("gr"));
               }catch (FileDoesNotExistException e){
                   System.out.println("Can not find the specified gene reference file, please check again.");
                   System.exit(-1);
               }
               GeneSet.readGeneInfo(saicPara.getGeneRefDir());
           }                      
           if (!cmd.hasOption("ta")){
               System.out.println("User didn't specify the threshold for copy number amplification detection, the software will use the default value" );
           }
           else{ 
               double ta = 0.0;
               try {
                    ta = Double.parseDouble(cmd.getOptionValue("ta"));
                    if (ta<0)
                    {
                        System.out.println("Invalid value for copy number deletion threshold. Please choose a positive float number");                     
                        System.exit(-1);
                    }
                }catch (NumberFormatException e)
                {
                    System.out.println("Invalid value for copy number deletion threshold. Please choose a positive float number");
                    System.exit(-1);
                }
               saicPara.setTa(ta);
           }
           if (!cmd.hasOption("td")){
               System.out.println("User didn't specify the threshold for copy number deletion detection, the software will use the default value");
           }
           else{ 
               double td = -0.0;
               try {
                    td = Double.parseDouble(cmd.getOptionValue("td"));
                    td = 0-td;
                    if (td>0)
                    {
                        System.out.println("Invalid value for copy number deletion threshold. Please choose a negative float number");                     
                        System.exit(-1);
                    }
                }catch (NumberFormatException e)
                {
                    System.out.println("Invalid value for copy number deletion threshold. Please choose a negative float number");
                    System.exit(-1);
                }
               saicPara.setTd(td);
           }
           
           if (!cmd.hasOption("np")){
               System.out.println("User didn't specify the number of permutation, the software will use the default value" );
           }
           else{ 
               int np = 0;
               try {
                    np = Integer.parseInt(cmd.getOptionValue("np"));
                    if (np<0)
                    {
                        System.out.println("Invalid value for number of permutation. Please choose in a positive integer." );                     
                        System.exit(-1);
                    }
                }catch (NumberFormatException e)
                {
                    System.out.println("Invalid value for number of permutation. Please choose in a positive integer.");
                    System.exit(-1);
                }
               saicPara.setNp(np);
           }
           
           if (!cmd.hasOption("tr")){
               System.out.println("User didn't specify the threshold for correltion coefficient between consecutive probes, the software will use the default value");
           }
           else{ 
               double tr = 0.0;
               try {
                    tr = Double.parseDouble(cmd.getOptionValue("tr"));
                    if (tr<0 || tr>1)
                    {
                        System.out.println("Invalid value for correlation coefficient threshold. Please type in a float number between 0 and 1 ");                     
                        System.exit(-1);
                    }
                }catch (NumberFormatException e)
                {
                    System.out.println("Invalid value for correlation coefficient threshold. Please type in a float number between 0 and 1 " );                     
                    System.exit(-1);
                }
               saicPara.setTci(tr);
           }
           if (!cmd.hasOption("ip")){
               System.out.println("User didn't specify whether to use SCA-exclusive permutation, the program will use this scheme as default");
           }
           else{
               if (cmd.getOptionValue("ip").equals("0"))
                   saicPara.setIterFlag(false);
               else
                   saicPara.setIterFlag(true);
           }
           if (!cmd.hasOption("qs")){
               System.out.println("User didn't specify whether to use quick SCAs detection, the program will adopt this as default");
           }
           else{
               if (cmd.getOptionValue("qs").equals("0"))
                   saicPara.setQuickSAIC(false);
               else
                   saicPara.setQuickSAIC(true);
           }
           if (saicPara.getAnalysisMode() == 2){
                String infiles[] = new String[1];
                String outfiles[] = new String[1];
                infiles[0] = saicPara.getInfileDir()+  "GenomeWide";
                outfiles[0] = saicPara.getOutfileDir() + "GenomeWideResult";
                Long startTime = System.currentTimeMillis();
                DetectSCAs detectSCAs = new DetectSCAs(saicPara, infiles, outfiles);                
                ForkJoinPool detectSCAForkJoinPool = new ForkJoinPool();                                       
                detectSCAForkJoinPool.invoke(detectSCAs); 
                System.out.println("The analyzing time for the whole genome detection is " + (System.currentTimeMillis()-startTime)/1000 + " seconds.");
                detectSCAs = null;
                detectSCAForkJoinPool = null;
           }
           else if (saicPara.getAnalysisMode() == 1){
                String infiles[] = new String[1];
                String outfiles[] = new String[1];
               for (int chr =22; chr>0; chr--) 
                {
                    infiles[0] = saicPara.getInfileDir()+ "Chr" + chr;
                    outfiles[0] = saicPara.getOutfileDir()+ "ResultChr" + chr;                                                         
                    DetectSCAs detectSCAs = new DetectSCAs(saicPara, infiles, outfiles);                                   
                    ForkJoinPool detectSCAForkJoinPool = new ForkJoinPool();                    
                    Long startTime = System.currentTimeMillis();    
                    System.out.println("Start detecting significant copy number amplifications and deletions parallelely");
                    detectSCAForkJoinPool.invoke(detectSCAs);                                        
                    System.out.println("The analyzing time for chromosome" + chr + " is " + (System.currentTimeMillis()-startTime)/1000 + " seconds.");
                    detectSCAs = null;
                    detectSCAForkJoinPool = null;
                }       
           }
           
           // Plot SAIC result figures

            double[] cn;
            String[] chrID;
            int[] loc;
            String fileName;
            double[] meanVal;
            boolean[] isSignificant;

            final String DATA_PATH = saicPara.getOutfileDir() + File.separator;
            final String AMP_FILE_SUFFIX = "_ampAvgCNA";
            final String DEL_FILE_SUFFIX = "_delAvgCNA";

            final String FIG_PATH = DATA_PATH + "figures" + File.separator;

            File figureDir = new File(FIG_PATH);

            // if the directory does not exist, create it
            if (!figureDir.exists()) {
                figureDir.mkdir();  
            }

            if ( saicPara.getAnalysisMode() == 2){
                final String FILE_PREFIX_GENOME = "ResultGenomeWide";
                String ampFile = DATA_PATH + FILE_PREFIX_GENOME +  AMP_FILE_SUFFIX;
                String delFile = DATA_PATH + FILE_PREFIX_GENOME + DEL_FILE_SUFFIX;
                String figureFileName = FIG_PATH + "AISAIC_genome.png";
                String genome = "genome";
                AisaicPlot p = createAisaicPlotObjectFromDataFiles(genome,
                        ampFile, delFile);
                p.plot();
                p.saveChartAsPNG(figureFileName);
            }

            else if(saicPara.getAnalysisMode() == 1){
                String[] chrArray = {"1", "2", "3", "4", "5", "6", "7", "8", 
                "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", 
                "19", "20", "21", "22"};
                for (String chr: chrArray) {
                    final String FILE_PREFIX_CHR = "ResultChr";
                    String ampFile = DATA_PATH + FILE_PREFIX_CHR + chr + AMP_FILE_SUFFIX;
                    String delFile = DATA_PATH + FILE_PREFIX_CHR + chr + DEL_FILE_SUFFIX;
                    String figureFileName = FIG_PATH + "AISAIC_chr_" + chr + ".png";
                    AisaicPlot p = createAisaicPlotObjectFromDataFiles(chr,
                            ampFile, delFile);
                    p.plot();
                    p.saveChartAsPNG(figureFileName);
                }
            }

       } 
        
    }
    
}
