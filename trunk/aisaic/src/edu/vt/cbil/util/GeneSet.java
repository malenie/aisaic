/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.vt.cbil.util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author bell
 */
public class GeneSet {
    
   private static ArrayList<Integer> geneChr ;
   
   private static ArrayList<Integer> geneStart;
   
   private static ArrayList<Integer> geneEnd;
   
   private static ArrayList<String> gene;
   
   public GeneSet()
   {
        geneChr = new ArrayList<Integer>(30000);
        
        geneStart = new ArrayList<Integer>(30000);
        
        geneEnd = new ArrayList<Integer>(30000);
        
        gene = new ArrayList<String>(30000);
   }
    
    
    
     public static void readGeneInfo (String geneRefDir)
    {
        geneChr = new ArrayList<Integer>(30000);
        
        geneStart = new ArrayList<Integer>(30000);
        
        geneEnd = new ArrayList<Integer>(30000);
        
        gene = new ArrayList<String>(30000);
        
        BufferedReader geneRefReader ;
        
        try{
            
            geneRefReader = new BufferedReader(new FileReader(geneRefDir));
            String temp;
            String[] splitStr;
            
            temp = geneRefReader.readLine();  //The first line is column names
            while((temp = geneRefReader.readLine())!=null)
            {
                splitStr = temp.split("\t");
                try{
                    geneChr.add(Integer.parseInt(splitStr[0]));
                }catch(NumberFormatException nfException)
                {
                    continue;                // continue or break depends on the input reference gene file
                }
                geneStart.add(Integer.parseInt(splitStr[1]));
                geneEnd.add(Integer.parseInt(splitStr[2]));
                gene.add(splitStr[4]);
            }

        }catch(IOException ioException)
        {
            System.out.println("Error reading file!");
        }
    }

    public static ArrayList<Integer> getGeneChr() {
        return geneChr;
    }

    public static void setGeneChr(ArrayList<Integer> geneChr) {
        GeneSet.geneChr = geneChr;
    }

    public static ArrayList<Integer> getGeneStart() {
        return geneStart;
    }

    public static void setGeneStart(ArrayList<Integer> geneStart) {
        GeneSet.geneStart = geneStart;
    }

    public static ArrayList<Integer> getGeneEnd() {
        return geneEnd;
    }

    public static void setGeneEnd(ArrayList<Integer> geneEnd) {
        GeneSet.geneEnd = geneEnd;
    }

    public static ArrayList<String> getGene() {
        return gene;
    }

    public static void setGene(ArrayList<String> gene) {
        GeneSet.gene = gene;
    }
     
     
    
}
