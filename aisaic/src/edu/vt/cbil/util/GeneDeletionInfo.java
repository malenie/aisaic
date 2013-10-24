/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.vt.cbil.util;

/**
 *
 * @author bell
 */
public class GeneDeletionInfo {
    
    private static int numHomoDeletion = 0;
    
    private static int numHemiDeletion = 0;
    
    private static GeneDeletionInfo instanceGeneDeletionInfo= null;
    
    private GeneDeletionInfo(){
       
    }
    
    public static synchronized GeneDeletionInfo getInstance()
    {        
        if (instanceGeneDeletionInfo==null)
        {
            instanceGeneDeletionInfo = new GeneDeletionInfo();
        }
        return instanceGeneDeletionInfo;
    }
    
    public synchronized void addHomoDeletion()
    {
        numHomoDeletion = numHomoDeletion + 1;
    }
    
    public synchronized void addHemiDeletion()
    {
        numHemiDeletion = numHemiDeletion + 1;
    }
    
    public int getNumHomoDeletion()
    {
        return numHomoDeletion;
    }
    
    public int getNumHemiDeletion()
    {
        return numHemiDeletion;
    }

    public void reset()
    {
         numHomoDeletion = 0;
         numHemiDeletion = 0;
    }
}
