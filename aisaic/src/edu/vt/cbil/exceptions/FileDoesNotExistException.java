/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.vt.cbil.exceptions;

/**
 *
 * @author houxuchu
 */
public class FileDoesNotExistException extends Exception{
    
    public FileDoesNotExistException(String errorMsg){
        super(errorMsg);
    }
    
}
