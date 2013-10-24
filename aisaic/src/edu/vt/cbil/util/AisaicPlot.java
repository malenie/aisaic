import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;
import java.io.File;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Shape;
import java.awt.BasicStroke;
import java.awt.Desktop;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;


import org.jfree.chart.JFreeChart;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.ChartFactory;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.data.xy.XYSeries;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYDotRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.StandardChartTheme;
import org.jfree.ui.RectangleInsets;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.util.ShapeUtilities;

import org.jfree.chart.annotations.XYLineAnnotation;



public class AisaicPlot {
    private XYSeries copyNumber;
    private int cnLength;
    private String chr;
    private String figureFileName;
    private JFreeChart chart;

    public AisaicPlot(double[] cn, int[] loc, 
		   String chr, String fileName){
	this.copyNumber = new XYSeries("Copy number");
	this.cnLength = cn.length;
	this.chr = chr;
	this.figureFileName = fileName;
	for (int i = 0; i < cn.length; i ++) {
	    this.copyNumber.add(loc[i], cn[i]);
	}
	this.chart = null;
    }

    public void plot() {
	XYSeriesCollection dataset = new XYSeriesCollection();
	dataset.addSeries(copyNumber);
	this.chart = ChartFactory.createScatterPlot
	    (
	     "Chr " + chr,
	     "Chr " + chr,
	     "copy number",
	     dataset,
	     PlotOrientation.VERTICAL, // Plot Orientation
	     true,
	     true,
	     false
	     );

	StandardChartTheme theme = new StandardChartTheme("JFree");
	theme.setPlotBackgroundPaint(Color.WHITE);
	theme.setRangeGridlinePaint(Color.GRAY);
	theme.setDomainGridlinePaint(Color.GRAY);
	theme.apply(chart);

	XYPlot xyPlot = (XYPlot) chart.getPlot();

	
	// Shape cross = ShapeUtilities.createDiagonalCross(3, 1);
        // renderer.setBaseShape(cross);
        // renderer.setBasePaint(Color.GRAY);
        // //changing the Renderer to XYDotRenderer
        // //xyPlot.setRenderer(new XYDotRenderer());
        XYDotRenderer xydotrenderer = new XYDotRenderer();
        xyPlot.setRenderer(xydotrenderer);
	//	xydotrenderer.setBasePaint(Color.GRAY, true);

	XYItemRenderer renderer = xyPlot.getRenderer();

        // xydotrenderer.setSeriesShape(0, cross);

        renderer.setSeriesPaint(0, Color.GRAY);
        // NumberAxis domain = (NumberAxis) xyPlot.getDomainAxis();
        // domain.setRange(0.00, 1e5);
        // domain.setTickUnit(new NumberTickUnit(1e4));
        // domain.setVerticalTickLabels(true);
        NumberAxis range = (NumberAxis) xyPlot.getRangeAxis();
        range.setRange(-1.0, 5.0);
        range.setTickUnit(new NumberTickUnit(0.5));

	chart.setBackgroundPaint(Color.WHITE);


    }

    public void plotSegMean(List<Integer> start, List<Integer> end, 
			    List<Double> segMean) {
	XYLineAnnotation line = null;
	XYPlot xyPlot = (XYPlot) chart.getPlot();
	for (int i = 0; i < start.size(); i ++) {
	    Color c;
	    if (segMean.get(i) < 1.98) {
		c = Color.BLUE;
	    } else if (segMean.get(i) > 2.02) {
		c = Color.RED;
	    } else {
		c = Color.GREEN;
	    }
	    line = new XYLineAnnotation(start.get(i) + 1,
					segMean.get(i), 
					end.get(i) + 1, 
					segMean.get(i),
					new BasicStroke(3), c);
	    xyPlot.addAnnotation(line);
	}
        
	
    }

    public void saveChartAsPNG() {
	try {
	    ChartUtilities.saveChartAsPNG(new File(figureFileName), chart, 2500, 500);
	} catch (Exception e) {
	    System.err.println("Problem occurred creating chart.");
	    System.err.println(e.getMessage());
	    e.printStackTrace();
	}
    }
    
    public static void main(String[] args) {
	
	double[] cn = null;
	String chr = "1";
	String[] chrID = null;
	int[] loc = null;
	String fileName = null;
	double[] meanVal = null;

	BufferedReader reader = null;
	try {
	    String line;
	    reader = new BufferedReader(new FileReader("data.csv"));
	    int numLine = 0;
	    while ((line = reader.readLine()) != null) {
		    numLine ++;
	    }

	    cn = new double[numLine];
	    loc = new int[numLine];
	    meanVal = new double[numLine];
	    chrID = new String[numLine];

	    reader = new BufferedReader(new FileReader("data.csv"));
	    for (int i = 0; i < numLine; i ++ ) {
		line = reader.readLine();
		String[] elements = line.split(",");
		cn[i] = Double.parseDouble(elements[3]);
		loc[i] = i + 1;
		meanVal[i] = Double.parseDouble(elements[4]);
		chrID[i] = elements[1];
	    }
 
	} catch (IOException e) {
	    e.printStackTrace();
	} finally {
	    try {
		if (reader != null) {
		    reader.close();
		}
	    } catch (IOException ex) {
		ex.printStackTrace();
	    }
	}


	String[] chrList = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
			    "11", "12", "13", "14", "15", "16", "17", "18", "19", 
			    "20", "21", "22"};

	double[] chr_cn;
	int[] chr_loc;
	double[] chr_meanVal;
	AisaicPlot p;
	CNSeg cnSeg;
	for (int i = 0; i < chrList.length; i ++ ) {
	    int s = 0;
	    int e = 0;
	    int j = 0;
	    while (!chrID[j].equals(chrList[i])) {
		j ++;
	    }
	    s = j;
	    while (j < chrID.length && chrID[j].equals(chrList[i])) {
		j ++;
	    }
	    e = j - 1;
	    chr_cn = Arrays.copyOfRange(cn, s, e + 1);
	    chr_loc = Arrays.copyOfRange(loc, s, e + 1);
	    chr_meanVal = Arrays.copyOfRange(meanVal, s, e + 1);
	    chr = chrList[i];
	    fileName = "figures/CNA_chr_" + chr + ".png";
	    p = new AisaicPlot(chr_cn, chr_loc, chr, fileName);
	    System.out.println("Plotting chromosome " + chr);
	    p.plot();
	    cnSeg = new CNSeg();
	    cnSeg.meanValToSeg(chr_meanVal, chr_loc);
	    p.plotSegMean(cnSeg.start, cnSeg.end, cnSeg.segMean);
	    p.saveChartAsPNG();
	    
	}
	try {
	    chr = "1";
	    Desktop.getDesktop().open(new File("figures/CNA_chr_" + chr + ".png"));
	} catch (Exception e) {
	    System.err.println(e.getMessage());
	}
    }
}
