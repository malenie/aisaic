package edu.vt.cbil.bacom;

import java.util.*;
import static java.lang.Math.*;

public class Filters {

    double[] x;

    public Filters(double[] signal) {
    	x = signal;
    }

    public void medianFilter(int windowSize) {
	int n = windowSize;
	int l = x.length;
	double[] window = new double[n];
	double[] xFiltered = new double[l];
	
	for (int i = 0; i < l; i ++) {
	    for (int j = (int) floor((1.0 - n) / 2.0); j <= floor((n - 1.0) / 2.0); j++) {

		int k = j - (int) floor((1.0 - n) / 2.0);

		if (i + j < 0) {
		    window[k] = x[0];
		} else if (i + j >= l) {
		    window[k] = x[l - 1];
		} else {
		    window[k] = x[i + j];
		}
	    }
	    xFiltered[i] = median(window);
	}

	System.arraycopy(xFiltered, 0, x, 0, l);

    }

    double median(double[] y) {
	
	Arrays.sort(y);
	if ( y.length % 2 == 1 )  {
	    return y[(y.length - 1) / 2];
	} else {
	    return (y[y.length / 2 - 1] + y[y.length / 2]) / 2;
	}

    }

}