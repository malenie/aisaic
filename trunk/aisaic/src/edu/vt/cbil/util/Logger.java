/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.vt.cbil.util;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import javax.xml.crypto.Data;

/**
 *
 * @author bella
 */
public class Logger {
    
    private static PrintWriter logger;

    public static void setLogger(String outFile) throws IOException {
        
        try{
            logger = new PrintWriter(new FileWriter(outFile, true), true);
        }
        catch (IOException e){
            System.out.println("Error create the new FileWriter");
            throw new IOException();
        }
    }
    
    public static void logging(String string){
        logger.write(getTime() + " INFO: " + string + System.getProperty("line.separator") + System.getProperty("line.separator"));  
        logger.flush();
    }
    
    public static String getTime(){
        DateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        Date today = Calendar.getInstance().getTime();
        return df.format(today);                
    }
    
    public static void close(){
        logger.close();
    }

}
