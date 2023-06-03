/**
 * MIT License

Copyright (c) 2015 Undral Erdnetsogt, Thorsten Wagner

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

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
import ij.process.ImageProcessor;

import java.awt.Color;
import java.awt.Frame;
import java.awt.Polygon;
import java.util.ArrayList;
import java.util.Vector;

import de.biomedical_imaging.ij.shapeSmoothingSlow.ComplexNumber;
import de.biomedical_imaging.ij.shapeSmoothingSlow.MyUsefulMethods;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;
import ij.plugin.RoiEnlarger;
import ij.plugin.*;

/**
 * 
 * @author Undral Erdenetsogt
 *
 */
public class ShapeSmoothingUtil {	

	private boolean drawContours = false;
	private boolean fillObjects = true;
	private boolean blackBackground = false;
	private int minimumFD;
	
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
			coef = fourierEngine(blob.getOuterContour(), thresholdValue, thresholdIsPercentual, ip,output,false);
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
			
			ArrayList<Polygon> inner = blob.getInnerContours();
			for (Polygon polygon : inner) {
				coef = fourierEngine(polygon, thresholdValue, thresholdIsPercentual, ip,output,true);
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
		
	}
	
	public void setMinimumNumberOfFD(int minimum){
		minimumFD = minimum;
	}
	
	public void setDrawContours(boolean b){
		drawContours = b;
	}
	
	public void setFillObjects(boolean b) {
		fillObjects = b;
	}
	
	public void setBlackBackground(boolean b){
		blackBackground = b;
	}
	
	
	private double[] calc_fourier_coefficients(Polygon equiCont) {
		
		double[] contourPoints = new double[2 * equiCont.npoints];
		
		// Konturpunkte in die "richtige" Datenstrukturübertragen
		/*
		 * a[2*k] = Re[k], 
		 * a[2*k+1] = Im[k], 0<=k<n
		 */
		int j = 0;
		for(int i = 0; i < equiCont.npoints; i++) {
			contourPoints[j] = equiCont.xpoints[i];
			contourPoints[j+1] = equiCont.ypoints[i];
			j=j+2;
		}
		DoubleFFT_1D ft = new DoubleFFT_1D(equiCont.npoints);

				
		// Fourier-Hintransformation
		ft.complexForward(contourPoints);	
		
		return contourPoints;
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
	private double[] fourierEngine(Polygon contourPolygon, double thresholdValue, boolean thresholdIsPercentual, ImageProcessor ip, boolean output, boolean isInner) {
		
		Polygon equiCont = toEquidistantPolygon(contourPolygon);
		int numOfContourPoints = equiCont.npoints;

		if (thresholdIsPercentual) {
			thresholdValue = (thresholdValue * numOfContourPoints) / 100;
			thresholdValue = thresholdValue<minimumFD?minimumFD:thresholdValue;
		}
		// thresholdValue darf nicht die Anzahl der Konturpunkte übersteigen!
		if (thresholdValue > numOfContourPoints) {
			thresholdValue = numOfContourPoints;
		}
			
	
		double[] all_coefficients = calc_fourier_coefficients(equiCont);
		
		
		// Filterung
		double[] filtered_coefficients = new double[(int)thresholdValue*2+1];

		for(int i = 0; i < (int)thresholdValue*2; i=i+2){
			filtered_coefficients[i] = all_coefficients[i];
			filtered_coefficients[i+1] = all_coefficients[i+1];
		}
		
		
	
		int loopFrom = (int)thresholdValue;
		if (loopFrom % 2 != 0) {
			loopFrom = loopFrom + 1;
		}
		int loopUntil = all_coefficients.length - (int)thresholdValue;
		if (loopUntil % 2 != 0) {
			loopUntil = loopUntil + 1;
		}
		
		for (int i = loopFrom; i < loopUntil; i++) {
			all_coefficients[i] = 0;
		}
		
			
		// Rücktransformation
		DoubleFFT_1D ft = new DoubleFFT_1D(equiCont.npoints);
		ft.complexInverse(all_coefficients, true);
		int[] xpoints = new int[all_coefficients.length/2];
		int[] ypoints = new int[all_coefficients.length/2];
		
		int j=0;
		for(int i = 0; i < all_coefficients.length; i=i+2) {
			xpoints[j] = (int) Math.round(all_coefficients[i]);
			ypoints[j] = (int) Math.round(all_coefficients[i+1]);
			j++;
		}
		
		if(blackBackground){
			ip.setValue(255);
			if(isInner){
				ip.setValue(0);
			}
			
		}else
		{
			ip.setValue(0);
			if(isInner){
				ip.setValue(255);
			}
		}
		
	
		Polygon poly = new Polygon(xpoints, ypoints, j);
		// Zeichnen
		Polygon pl = null;
		int enlarge=0;
		int roiType = Roi.TRACED_ROI;
		if(drawContours == true & fillObjects == false) {
			ip.drawPolygon(poly);
		}
		else if(drawContours == true & fillObjects == true){
			ip.fillPolygon(poly);
			ip.setValue(255);
			if(blackBackground) {
				ip.setValue(0);
			}
			try {
				// This try catch is necessary as it otherwise crashes during preview in some cases.
				Roi roi = new PolygonRoi(poly,roiType);
				Roi newRoi = RoiEnlarger.enlarge(roi, enlarge);
				pl = newRoi.getPolygon();
				ip.drawPolygon(pl);
			} catch(NullPointerException e) {
				ip.drawPolygon(poly);
			}
		}
		else if(drawContours == false & fillObjects == true) {
			ip.drawPolygon(poly);
			ip.fillPolygon(poly);
		}
		
		
		if(output){
			Frame frame = WindowManager.getFrame("ROI Manager");
			if (frame == null)
				IJ.run("ROI Manager...");
			frame = WindowManager.getFrame("ROI Manager");
			RoiManager roiManager = (RoiManager) frame;
			Roi roi = new PolygonRoi(poly,roiType);
			roiManager.add(IJ.getImage(), roi,objectCounter);
			objectCounter++;
		}
		
		return filtered_coefficients;
	}
	
	public void drawPolygon(ImageProcessor ip, Polygon p) {
		
		ip.moveTo(p.xpoints[0], p.ypoints[0]);
		for (int i=0; i<p.npoints; i++)
			ip.lineTo(p.xpoints[i], p.ypoints[i]);
		ip.lineTo(p.xpoints[0], p.ypoints[0]);
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
