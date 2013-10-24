/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.vt.cbil.exceptions;

/**
 *
 * @author bell
 */
public class TooSmallSampleSizeException extends Exception {
    
    public TooSmallSampleSizeException(String errorMessage){
        super(errorMessage);
    }
    
}
