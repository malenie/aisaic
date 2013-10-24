/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.vt.cbil.saic;

import edu.vt.cbil.util.SaicParameters;
import java.util.concurrent.ForkJoinPool;

/**
 *
 * @author houxuchu
 * In order to call SAIC analysis, the user only need to instantiate a SaicPlugin object
 * and call compute() method
 */
public class SaicPlugin {
    
    public void compute(SaicParameters para, int type, double input[][]){        
        SAIC saic = new SAIC (para, type, input);
        ForkJoinPool saicForkJoinPool = new ForkJoinPool();
        saicForkJoinPool.invoke(saic);   
    }   
}
