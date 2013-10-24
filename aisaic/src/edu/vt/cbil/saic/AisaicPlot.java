/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.vt.cbil.saic;

import java.util.List;
import java.io.File;
import java.awt.Color;
import java.awt.Shape;
import java.awt.BasicStroke;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;


import org.jfree.chart.JFreeChart;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.ChartFactory;
import org.jfree.data.xy.XYSeries;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.chart.StandardChartTheme;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.util.ShapeUtilities;

import org.jfree.chart.annotations.XYLineAnnotation;


/**
 *
 * @author bai
 */
public class AisaicPlot {
    private XYSeries seriesAmpSig;
    private XYSeries seriesAmpInsig;
    private XYSeries seriesDelSig;
    private XYSeries seriesDelInsig;
    private int cnLength;
    private String chr;
    private String figureFileName;
    private JFreeChart chart;

    public AisaicPlot(double[] cnAmp, double[] cnDel, int[] loc, 
		      boolean[] isSignificantAmp, boolean[] isSignificantDel,
		      String chr) {
	this.seriesAmpSig = new XYSeries("significant amplification");
	this.seriesDelSig = new XYSeries("significant deletion");
	this.seriesAmpInsig = new XYSeries("insignificant amplification");
	this.seriesDelInsig = new XYSeries("insignificant deletion");
	this.cnLength = loc.length;
	this.chr = chr;

	this.seriesAmpInsig.add(0, 0);
	this.seriesDelInsig.add(0, 0);

	for (int i = 0; i < cnAmp.length; i ++) {
	    if (isSignificantAmp[i]) {
		this.seriesAmpSig.add(loc[i], cnAmp[i]);
	    } else {
		this.seriesAmpInsig.add(loc[i], cnAmp[i]);
	    }
	}
	for (int i = 0; i < cnDel.length; i ++) {
	    if (isSignificantDel[i]) {
		this.seriesDelSig.add(loc[i], cnDel[i]);
	    } else {
		this.seriesDelInsig.add(loc[i], cnDel[i]);
	    }
	}
	this.chart = null;
    }

    public void plot() {
	XYSeriesCollection dataset = new XYSeriesCollection();
	dataset.addSeries(seriesAmpSig);
	dataset.addSeries(seriesAmpInsig);
	dataset.addSeries(seriesDelSig);
	dataset.addSeries(seriesDelInsig);
	this.chart = ChartFactory.createScatterPlot
	    (
	     "Chr " + chr,
	     "Chr " + chr,
	     "copy number (log2)",
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

	

	XYItemRenderer renderer = xyPlot.getRenderer();
	renderer.setSeriesShape(0, (Shape) ShapeUtilities.createDiamond(2));
	renderer.setSeriesShape(1, (Shape) ShapeUtilities.createDiamond(1));
	renderer.setSeriesShape(2, (Shape) ShapeUtilities.createDiamond(2));
	renderer.setSeriesShape(3, (Shape) ShapeUtilities.createDiamond(1));

        renderer.setSeriesPaint(0, new Color(215, 25, 28));
        renderer.setSeriesPaint(1, new Color(253, 174, 97));
        renderer.setSeriesPaint(2, new Color(44, 123, 182));
        renderer.setSeriesPaint(3, new Color(171, 217, 233));

        NumberAxis range = (NumberAxis) xyPlot.getRangeAxis();
        range.setRange(-2.0, 2.0);
        range.setTickUnit(new NumberTickUnit(0.5));

	chart.setBackgroundPaint(Color.WHITE);


    }

   
    public void saveChartAsPNG(String figureFileName) {
	this.figureFileName = figureFileName;
	try {
	    ChartUtilities.saveChartAsPNG(new File(figureFileName), chart, 2500, 500);
	} catch (Exception e) {
	    System.err.println("Problem occurred creating chart.");
	    System.err.println(e.getMessage());
	    e.printStackTrace();
	}
    }


    public static AisaicPlot createAisaicPlotObjectFromDataFiles(String chr,
            String ampFile, String delFile) {

	double[] cnAmp = null;
	double[] cnDel = null;
	int[] loc = null;
	boolean[] isSignificantAmp = null;
	boolean[] isSignificantDel = null;
        


	BufferedReader reader = null;
	try {
	    String line;
	    reader = new BufferedReader(new FileReader(ampFile));
	    int numLine = 0;
	    while ((line = reader.readLine()) != null) {
		    numLine ++;
	    }

	    cnAmp = new double[numLine];
	    loc = new int[numLine];
	    isSignificantAmp = new boolean[numLine];

	    reader = new BufferedReader(new FileReader(ampFile));
	    for (int i = 0; i < numLine; i ++ ) {
		line = reader.readLine();
		String[] elements = line.split("\t");
		cnAmp[i] = Double.parseDouble(elements[0]);
		loc[i] = i + 1;
		isSignificantAmp[i] = Integer.parseInt(elements[1]) == 0? false : true;
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


	try {
	    String line;
	    reader = new BufferedReader(new FileReader(delFile));
	    int numLine = 0;
	    while ((line = reader.readLine()) != null) {
		    numLine ++;
	    }

	    cnDel = new double[numLine];
	    isSignificantDel = new boolean[numLine];

	    reader = new BufferedReader(new FileReader(delFile));
	    for (int i = 0; i < numLine; i ++ ) {
		line = reader.readLine();
		String[] elements = line.split("\t");
		cnDel[i] = Double.parseDouble(elements[0]);
		isSignificantDel[i] = Integer.parseInt(elements[1]) == 0? false : true;
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

        return new AisaicPlot(cnAmp, cnDel, loc, 
                isSignificantAmp, isSignificantDel,
                chr);

    }

          
}
