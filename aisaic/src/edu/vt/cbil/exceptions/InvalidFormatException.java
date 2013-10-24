/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.vt.cbil.exceptions;

/**
 *
 * @author bell
 */
public class InvalidFormatException extends Exception{
    
    public InvalidFormatException(String errorMessage){
        super(errorMessage);
    }
}
