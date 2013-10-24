import java.util.ArrayList;

public class CNSeg {

    public ArrayList<Integer> start;
    public ArrayList<Integer> end;
    public ArrayList<Double> segMean;
	
    public CNSeg() {
	start = new ArrayList<Integer>() ;
	end = new ArrayList<Integer>();
	segMean = new ArrayList<Double>();
    }

    public void meanValToSeg(double[] meanVal) {
	int currentStart = 0;
	double currentMean = meanVal[0];
	start.add(currentStart);
	segMean.add(currentMean);
	for (int i = 1; i < meanVal.length; i ++) {
	    if (meanVal[i] < currentMean - 0.001 || meanVal[i] > currentMean + 0.001) {
		end.add(i - 1);
		currentStart = i;
		currentMean = meanVal[i];
		start.add(currentStart);
		segMean.add(currentMean);
	    }
	}
	end.add(meanVal.length - 1);
    }

    public void meanValToSeg(double[] meanVal, int[] loc) {
	int currentStart = loc[0];
	double currentMean = meanVal[0];
	start.add(currentStart);
	segMean.add(currentMean);
	for (int i = 1; i < meanVal.length; i ++) {
	    if (meanVal[i] < currentMean - 0.001 || meanVal[i] > currentMean + 0.001) {
		end.add(loc[i - 1]);
		currentStart = loc[i];
		currentMean = meanVal[i];
		start.add(currentStart);
		segMean.add(currentMean);
	    }
	}
	end.add(loc[meanVal.length - 1]);
    }

}