/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.vt.cbil.util;

import edu.vt.cbil.exceptions.InvalidFormatException;

/**
 *
 * @author Admin
 */
public class Gene {
    
    private String chr;
    private int start;
    private int end;

    public String getChr() {
        return chr;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }
    
    public void parseStrToGene(String input) throws InvalidFormatException
    {   
        String[] splitStr;
         
         try{
             splitStr = input.split(":");
         }
         catch (Exception e){
             throw new InvalidFormatException ("Invalid input");
         }
         
         if (splitStr.length != 2){
             throw new InvalidFormatException("Invalid input");
         }
         chr = splitStr[0]; 
         
         try{
             Integer.parseInt(chr);
         }
         catch (Exception e){
             throw new InvalidFormatException("Invalid input");
         }
         
         String loci = splitStr[1];
         
         try{
            splitStr = loci.split("-");
         }
         catch (Exception e){
             throw new InvalidFormatException ("Invalid input");
         }
                  
         if (splitStr.length != 2){
             throw new InvalidFormatException("Invalid input");
         }
         
         try{
             start = Integer.parseInt(splitStr[0]); 
         }
         catch (Exception e){
             throw new InvalidFormatException("Invalid input");
         }
         
         try{
            end = Integer.parseInt(splitStr[1]); 
         }
         catch (Exception e){
             throw new InvalidFormatException("Invalid input");
         }
         
        
    }
    
}
