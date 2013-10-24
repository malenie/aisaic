/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.vt.cbil.util;

import edu.vt.cbil.exceptions.EmptyStringException;
import edu.vt.cbil.exceptions.FileDoesNotExistException;
import edu.vt.cbil.exceptions.InvalidFormatException;
import edu.vt.cbil.exceptions.InvalidInputFileException;
import edu.vt.cbil.exceptions.SampleSizeDoesNotMatchException;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author houxuchu
 */
public class BacomParameters {
    
    private String analysisMode;
    private String inputNormalFile;
    private String inputTumorFile;
    private String inputCdfFile;
    private String annotationFile;
    private String resultFolder;
    private Gene   geneInterested;
    private ArrayList<String> sampleIds;
    private ArrayList<String> normalSamples;
    private ArrayList<String> tumorSamples;
    private double [] alpha;
    private int sampleSize;

    public BacomParameters(String analysisMode, String inputNormalFile, String inputTumorFile, String resultFolder, String geneInterested) 
                      throws InvalidFormatException,EmptyStringException{
        setAnalysisMode(analysisMode);
        setInputNormalFile(inputNormalFile);
        setInputTumorFile (inputTumorFile);
        setResultFolder (resultFolder);
        setGeneInterested(geneInterested);
        
    }

    public ArrayList<String> getSampleIds() {
        return sampleIds;
    }

    public ArrayList<String> getNormalSamples() {
        return normalSamples;
    }

    public ArrayList<String> getTumorSamples() {
        return tumorSamples;
    }

    public double[] getAlpha() {
        return alpha;
    }

    public int getSampleSize() {
        return sampleSize;
    }

    public void setNormalSamples() throws InvalidInputFileException, IOException, FileDoesNotExistException{
        if (analysisMode.equals("0")){
            int dotIndex = inputNormalFile.lastIndexOf(".");
            String fileType = inputNormalFile.substring(dotIndex);  
            if (".cel".equals(fileType.toLowerCase())){
                normalSamples = new ArrayList<String> (1);
                File tempFile = new File(inputNormalFile);
                if (!tempFile.exists()){
                    throw new FileDoesNotExistException("File does not exist, please check.");
                }
                normalSamples.add(inputNormalFile);  
                sampleSize = 1;
            }
            else {
                throw new InvalidInputFileException("Invalid input file!");
            }
        }
        else if (analysisMode.equals("1")){  
            normalSamples = new ArrayList<String> (500);
            BufferedReader normalSamplesBuffer = new BufferedReader(new FileReader(inputNormalFile));                                  
            String normalTemp;                    
            String normalFileDir = new File(inputNormalFile).getParent();
            int i=0;   
            File tempFile;
            while ((normalTemp = normalSamplesBuffer.readLine())!=null)
            {
                tempFile = new File(normalFileDir + File.separator + normalTemp);
                if (!tempFile.exists())
                    throw new FileDoesNotExistException("File does not exist, please check.");
                normalSamples.add(tempFile.getAbsolutePath());                                                
                i++;
            }   
            sampleSize = i;
            normalSamples.trimToSize();                
        }
         alpha = new double[sampleSize];  
    }

    public void setTumorSamples() throws SampleSizeDoesNotMatchException, InvalidInputFileException, IOException, FileDoesNotExistException{
        if (analysisMode.equals("0")){
            int dotIndex = inputTumorFile.lastIndexOf(".");
            String fileType = inputTumorFile.substring(dotIndex);
            if (".cel".equals(fileType.toLowerCase())){
                File tempFile = new File(inputTumorFile);  
                if (!tempFile.exists()){
                    throw new FileDoesNotExistException("File does not exist, please check.");
                }
                tumorSamples = new ArrayList<String> (1);
                tumorSamples.add(inputTumorFile);   
                sampleIds = new ArrayList<String> (1);                        
                sampleIds.add(tempFile.getName());  
            }
            else{
                throw new InvalidInputFileException("Invalid input file!");
            }
        }
        else if (analysisMode.equals("1")){  
            tumorSamples = new ArrayList<String> (500);
            sampleIds = new ArrayList<String> (500);
            try{
                BufferedReader tumorSampleBuffer = new BufferedReader(new FileReader(inputTumorFile));                                  
                String tumorTemp;                    
                String tumorFileDir = new File(inputTumorFile).getParent();  
                int i=0;
                File tempFile;
                while ((tumorTemp = tumorSampleBuffer.readLine())!=null)
                {
                    sampleIds.add(tumorTemp);
                    tempFile = new File(tumorFileDir + File.separator + tumorTemp);
                    if (!tempFile.exists()){
                        throw new FileDoesNotExistException("File does not exist, please check.");
                    }
                    tumorSamples.add(tempFile.getAbsolutePath());  
                    i++;
                }   
                if (sampleSize != i){
                    throw new SampleSizeDoesNotMatchException("Sample size does not MATCH! Please check your input");       
                }              
                tumorSamples.trimToSize();                              
            }
            catch (Exception e){
                throw e;
            }
        }
    }
    
    public void setGeneInterested(String geneString) throws InvalidFormatException, EmptyStringException{
        if (geneString == null || geneString.equals("")){
            throw new EmptyStringException("Empty string!");
        }
        try{
            geneInterested = new Gene();
            geneInterested.parseStrToGene(geneString);
        }
        catch (Exception e){
            throw new InvalidFormatException("Invalid input");
        }
    }
    
    public void setGeneInterestedNull(){
        geneInterested = null;
    }


    public BacomParameters() {
        analysisMode = "0";         
    }

    public String getAnalysisMode() {
        return analysisMode;
    }

    public String getInputNormalFile() {
        return inputNormalFile;
    }

    public String getInputTumorFile() {
        return inputTumorFile;
    }

    public String getResultFolder() {
        return resultFolder;
    }

    public Gene getGeneInterested() {
        return geneInterested;
    }

    public void setAnalysisMode(String analysisMode) {
        this.analysisMode = analysisMode;
    }

    public void setInputNormalFile(String inputNormalFile) {
        this.inputNormalFile = inputNormalFile;
    }

    public void setInputTumorFile(String inputTumorFile) {
        this.inputTumorFile = inputTumorFile;
    }

    public void setResultFolder(String resultFolder) {
        this.resultFolder = resultFolder;
    }

    public void setInputCdfFile(String inputCdfFile) throws FileDoesNotExistException{
        File cdfFile = new File(inputCdfFile);
        if(!cdfFile.exists())
            throw new FileDoesNotExistException("The specified CDF file does not exist, please check!");
        this.inputCdfFile = inputCdfFile;
    }

    public String getInputCdfFile() {
        return inputCdfFile;
    }

    public String getAnnotationFile() {
        return annotationFile;
    }

    public void setAnnotationFile(String annotationFile)throws FileDoesNotExistException {
        File annoFile = new File(annotationFile);
        if(!annoFile.exists()){
            throw new FileDoesNotExistException("The specified annotation file does not exist, please check!");
        }
        this.annotationFile = annotationFile;
    }
    
    
}
