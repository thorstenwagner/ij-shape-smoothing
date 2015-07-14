package de.biomedical_imaging.ij.shapeSmoothing;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.blob.Blob;
import ij.blob.ManyBlobs;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.measure.ResultsTable;
import ij.plugin.frame.RoiManager;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.awt.Frame;
import java.awt.Polygon;
import java.util.Vector;

import de.biomedical_imaging.ij.shapeSmoothingPlugin.Shape_Smoothing;
import de.biomedical_imaging.ij.shapeSmoothingSlow.ComplexNumber;
import de.biomedical_imaging.ij.shapeSmoothingSlow.MyUsefulMethods;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;

/**
 * 
 * @author Undral Erdenetsogt
 *
 */
public class ShapeSmoothingUtil {	

	private boolean onlyContours = false;
	private boolean blackBackground = false;
	/**
	 * Fourier-Hintransformation, Filterung der Fourierdeskriptoren (FD) und Fourier-Rücktransformation
	 * für alle Blobs auf einem Bild (referenziert mittels imp). 
	 * 
	 * @param imp {@link ImagePlus} Das Bild auf dem nach Blobs gesucht wird. Vorsicht, Inhalte werden überschrieben!
	 * @param thresholdValue Schwellenwert für die Filterung der FDs - es gibt an, wie viele FDs beibehalten werden (absolut oder prozentual)
	 * @param thresholdIsPercentual Gibt an, ob thresholdValue eine prozentuale Angabe ist, oder ob es bereits die Anzahl der zu behaltenden FDs enthält 
	 * @param output Ausgabe der Descriptoren in Result Table
	 */
	
	public void fourierFilter(ImageProcessor ip, double thresholdValue, boolean thresholdIsPercentual, boolean output) {
		// Konturpunkte erfassen	
		ImagePlus imp = new ImagePlus("", ip);
		ManyBlobs allBlobs = new ManyBlobs(imp); // Extended ArrayList
		if(blackBackground){
			allBlobs.setBackground(0);
		}else {
			allBlobs.setBackground(1);
			
		}
		allBlobs.findConnectedComponents(); // Formen erkennen	
		
		if(blackBackground){
			ip.setColor(Color.black);
		}else{
			ip.setColor(Color.white);
		}
		ip.fill();
		ResultsTable rt = new ResultsTable();
		rt.setPrecision(5);
		
		// Clean ROI Manager
		Frame frame = WindowManager.getFrame("ROI Manager");
		if (frame != null){
			RoiManager roiManager = (RoiManager) frame;
			roiManager.close();
		}
		int c = 1;
		for (Blob blob: allBlobs) {
			double[] coef;
			coef = fourierEngine(blob.getOuterContour(), thresholdValue, thresholdIsPercentual, ip,output);
			if(output && coef.length >= 4){
				double f1 = Math.sqrt(Math.pow(coef[2],2)+Math.pow(coef[3],2));
				rt.incrementCounter();
				rt.addValue("Blob Label", c);
				for(int i = 0; i< coef.length-1; i=i+2){
					
					rt.showRowNumbers(false);
					
					//rt.addValue("R", coef[i]);
				//	rt.addValue("I", coef[i+1]);
					rt.addValue("|F" +i/2+ "|", Math.sqrt(Math.pow(coef[i],2)+Math.pow(coef[i+1],2))/f1);
				}
				rt.show("Fourier Descriptors");
			}
			c++;
		}
		
	}
	

	
	public void setDrawOnlyContours(boolean b){
		onlyContours = b;
	}
	
	public void setBlackBackground(boolean b){
		blackBackground = b;
	}
	
	/**
	 * Fourier-Hintransformation, Filterung der Fourierdeskriptoren (FD) und Fourier-Rücktransformation
	 * für einen Blob (referenziert mittels contourPolygon). 
	 * 
	 * @param contourPolygon Kontur eines Blobs
	 * @param thresholdValue Schwellenwert für die Filterung der FDs (prozentual oder absolut): die ersten thresholdValue/2 und die letzten thresholdValue/2 FDs werden
	 * beibehalten und die Restlichen auf 0 gesetzt
	 * @param thresholdIsPercentual Gibt an, ob thresholdValue eine prozentuale Angabe ist, oder ob es bereits die Anzahl der zu behaltenden FDs enthält 
	 * @param ip {@link ImageProcessor} des Bildes, worauf das rücktransformierte Blob gezeichnet wird
	 */
	int objectCounter = 1;
	private double[] fourierEngine(Polygon contourPolygon, double thresholdValue, boolean thresholdIsPercentual, ImageProcessor ip, boolean output) {
		
		Polygon equiCont = toEquidistantPolygon(contourPolygon);
		int numOfContourPoints = equiCont.npoints;
		
		if (thresholdIsPercentual) {
			thresholdValue = (thresholdValue * numOfContourPoints) / 100;
			//
		} else {
			// thresholdValue darf nicht die Anzahl der Konturpunkte übersteigen!
			if (thresholdValue > numOfContourPoints) {
				thresholdValue = numOfContourPoints;
			}
		}
			
		// Konturpunkte in die "richtige" Datenstrukturübertragen

		/*
		 * a[2*k] = Re[k], 
		 * a[2*k+1] = Im[k], 0<=k<n
		 */
		double[] contourPoints = new double[2 * numOfContourPoints];
		
		int j = 0;
		for(int i = 0; i < numOfContourPoints; i++) {
			contourPoints[j] = equiCont.xpoints[i];
			contourPoints[j+1] = equiCont.ypoints[i];
			j=j+2;
		}
		DoubleFFT_1D ft = new DoubleFFT_1D(numOfContourPoints);

				
		// Fourier-Hintransformation
		ft.complexForward(contourPoints);		
		double[] coefficients = new double[(int)thresholdValue*2+1];
		
		for(int i = 0; i < (int)thresholdValue*2; i=i+2){
			coefficients[i] = contourPoints[i];
			coefficients[i+1] = contourPoints[i+1];
		}
		
		// Filterung	
	
		int loopFrom = (int)thresholdValue;
		if (loopFrom % 2 != 0) {
			loopFrom = loopFrom + 1;
		}
		int loopUntil = contourPoints.length - (int)thresholdValue;
		if (loopUntil % 2 != 0) {
			loopUntil = loopUntil + 1;
		}
		
		for (int i = loopFrom; i < loopUntil; i++) {
			contourPoints[i] = 0;
		}
		
			
		// Rücktransformation
		ft.complexInverse(contourPoints, true);
		int[] xpoints = new int[contourPoints.length/2];
		int[] ypoints = new int[contourPoints.length/2];
		
		j=0;
		for(int i = 0; i < contourPoints.length; i=i+2) {
			xpoints[j] = (int) Math.round(contourPoints[i]);
			ypoints[j] = (int) Math.round(contourPoints[i+1]);
			j++;
		}

		if(blackBackground){
			//IJ.setForegroundColor(255, 255, 255);
			ip.setValue(255);
			
		}else
		{
			ip.setValue(0);
			//IJ.setForegroundColor(0, 0, 0);
			
		}
		Polygon poly = new Polygon(xpoints, ypoints, j);
		// Zeichnen
		if(onlyContours){
			ip.drawPolygon(poly);
		}else{
			ip.drawPolygon(poly);
			ip.fillPolygon(poly);
		}
		
		if(output){
			Frame frame = WindowManager.getFrame("ROI Manager");
			if (frame == null)
				IJ.run("ROI Manager...");
			frame = WindowManager.getFrame("ROI Manager");
			RoiManager roiManager = (RoiManager) frame;
			Roi roi = new PolygonRoi(poly,Roi.TRACED_ROI);
			roiManager.add(IJ.getImage(), roi,objectCounter);
			objectCounter++;
		}
		
		return coefficients;
	}
	
	public Polygon toEquidistantPolygon(Polygon pol){
		Polygon equiPol = new Polygon();
	
		for(int i = 0; i < pol.npoints; i++){
			equiPol.addPoint(pol.xpoints[i], pol.ypoints[i]);
			int nextPixel = (i+1)==pol.npoints? 1 : i+1;
			int dx = pol.xpoints[i] - pol.xpoints[nextPixel];
			int dy = pol.ypoints[i] - pol.ypoints[nextPixel];
			if((dx == -1 && dy ==-1) || (dx == 1 && dy ==1)){
				equiPol.addPoint(pol.xpoints[nextPixel],  pol.ypoints[i] );
			}

			else if((dx == -1 && dy == 1) || (dx == 1 && dy == -1)){
				equiPol.addPoint(pol.xpoints[i],  pol.ypoints[nextPixel]);
			}
		}
		return equiPol;
	}
	
	// Es ist nur für Zeitmesszwecke da	
	@SuppressWarnings("unused")
	@Deprecated
	private void fourierEngineSlow (Polygon contourPolygon, int thresholdValue, boolean thresholdIsPercentual, ImageProcessor ip) {
	
		if (thresholdIsPercentual) {
			thresholdValue = (thresholdValue * contourPolygon.npoints) / 100;
		}
		
		// Konturpunkte in die "richtige" Datenstruktur übertragen
		Vector<ComplexNumber> contourPoints = new Vector<ComplexNumber>();
		for (int i = 0; i < contourPolygon.npoints; i++) {
			contourPoints.add(new ComplexNumber(contourPolygon.xpoints[i], contourPolygon.ypoints[i]));
		}
		
		
		// Fourier-Hintransformation
		Vector<ComplexNumber> FDs = MyUsefulMethods.dft(contourPoints);		
		
		// Filterung		
		int loopFrom = thresholdValue/2;
		if (thresholdValue % 2 != 0) {
			loopFrom = loopFrom + 1; // bei ungeradem threshold aufrunden;
		}
		int loopUntil = FDs.size() - thresholdValue/2; // wird 'von allein' abgerundet
		ComplexNumber complexZero = new ComplexNumber(0, 0);
		for (int i = loopFrom; i < loopUntil ; i++ ) {
			FDs.setElementAt(complexZero, i);
		}
		
		
		// Rücktransformation
		Vector<ComplexNumber> newPoints = MyUsefulMethods.dftInverse(FDs);
	
		
		// Zeichnen
		for (ComplexNumber newPoint: newPoints) {
			ip.putPixel((int) Math.round(newPoint.getRe()), (int) Math.round(newPoint.getIm()), 0);
		}

	}
		
	
}
