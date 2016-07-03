package de.biomedical_imaging.ij.shapeSmoothingPlugin;

import java.awt.AWTEvent;
import java.awt.Scrollbar;
import java.awt.TextField;
import java.awt.Choice;

import ij.IJ;
import ij.ImagePlus;
import ij.blob.ManyBlobs;
import ij.blob.Blob;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NewImage;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import de.biomedical_imaging.ij.shapeSmoothing.ShapeSmoothingUtil;

/**
 * Das ImageJ Plugin macht eine Konturglättung der Formen auf einem Bild. Das Plugin ist nach dem Pipeline-Prinzip entwickelt.
 * 
 * @author Undral Erdenetsogt
 * 
 * V1.1: Adding preview option (Thorsten Wagner)
 *
 */
public class Shape_Smoothing implements ExtendedPlugInFilter, DialogListener {

	private ShapeSmoothingUtil shapeSmoothingUtil;
	private double thresholdValuePercentual;
	private double thresholdValueAbsolute;
	private Choice modusChoice ;
	private boolean doAbsoluteThreshold;
	private boolean drawOnlyContours;
	private boolean blackBackground;
	private boolean doOutputDescriptors = false;
	String[] absRelChoices = {"Relative proportion of FDs","Absolute number of FDs"};
	int maxNumOfFDs ;
	int minNumOfFDs;
	int numberOfBlobs;
	public static boolean invertedLut;
	boolean previewing = true;
	@Override
	public int setup(String arg, ImagePlus imp) {
		
		//imp = ensureCorrectLUT(imp);
		shapeSmoothingUtil = new ShapeSmoothingUtil();
		invertedLut = imp.isInvertedLut();
		if (imp == null || imp.getType() != ImagePlus.GRAY8) {
			IJ.error("Only 8-Bit Grayscale Imags are supported");
			return DONE;
		}
		//YesNoCancelDialog diag = new YesNoCancelDialog(ImageWindow.getFrames()[0], "Background", "Black background?");
		ImageStatistics stats = imp.getStatistics();

		blackBackground = (stats.histogram[0]>stats.histogram[255]); //diag.yesPressed();
		//if(invertedLut) blackBackground = !blackBackground;

		
		ManyBlobs allBlobs = new ManyBlobs(imp);
		if(blackBackground){
			allBlobs.setBackground(0);
		}else
		{
			allBlobs.setBackground(1);
		}
		allBlobs.findConnectedComponents();
		numberOfBlobs = allBlobs.size();
		minNumOfFDs = Integer.MAX_VALUE;
		int tempMaxNumOfFDs = 0;
		
		for (Blob blob : allBlobs) {
			
			int numOfFDs = shapeSmoothingUtil.toEquidistantPolygon(blob.getOuterContour()).npoints;
			if (numOfFDs < minNumOfFDs) {
				minNumOfFDs = numOfFDs;
			}
			
			if (numOfFDs > tempMaxNumOfFDs) {
				tempMaxNumOfFDs = numOfFDs;
			}
		}
		
		maxNumOfFDs = tempMaxNumOfFDs;
		
		
		return DOES_8G;
	}
	
	@Override
	public void run(ImageProcessor ip) {		
		doFourierFilter(ip);
	}
	
	@Override
	public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
		IJ.showStatus("Fourier Hin- & Rücktransformation");
		
		GenericDialog gd = new GenericDialog(command + "...");
		gd.setOKLabel("Run");
		gd.setCancelLabel("Cancel");
		// gd.addMessage("There are " + numberOfBlobs + " objects (or contours) with " + minNumOfFDs + " to " + maxNumOfFDs + " Fourier Descriptors (FDs).");		
		gd.addChoice("Keep (for each blob):", absRelChoices, absRelChoices[0]);
		gd.addSlider("Relative proportion FDs (%)", 0, 100, 2);
		gd.addSlider("Absolute number FDs", 1, maxNumOfFDs, 2);
		//dialogItemChanged(gd, null);
		gd.addCheckbox("Draw only contours", false);
		gd.addCheckbox("Output Descriptors", false);

		gd.addCheckbox("Black Background", blackBackground);
		Scrollbar absScroll = (Scrollbar) gd.getSliders().get(1);
		TextField absTextField = (TextField) gd.getNumericFields().get(1);
		absScroll.setEnabled(false);
		absTextField.setEnabled(false);
		
		gd.addDialogListener(this);
		gd.addPreviewCheckbox(pfr);
		gd.showDialog();
		if (gd.wasCanceled()) {
			return DONE;
		}
		
		thresholdValuePercentual = gd.getNextNumber();
		thresholdValueAbsolute = (int) gd.getNextNumber();
		String choice = gd.getNextChoice();
		doAbsoluteThreshold = (choice.equals(absRelChoices[1]));
		drawOnlyContours = gd.getNextBoolean();
		doOutputDescriptors = gd.getNextBoolean();
		blackBackground = gd.getNextBoolean();
		previewing = false;
		return IJ.setupDialog(imp, DOES_8G);
	}
	
	int percent = 0;
	int absolute = 0;
	@Override
	public boolean dialogItemChanged(GenericDialog geDi, AWTEvent e) {
		thresholdValuePercentual = geDi.getNextNumber();
		thresholdValueAbsolute = (int) geDi.getNextNumber();
		modusChoice = (Choice) geDi.getChoices().get(0);
		doAbsoluteThreshold = (modusChoice.getSelectedItem()== absRelChoices[1]);
		drawOnlyContours = geDi.getNextBoolean();
		doOutputDescriptors = geDi.getNextBoolean();
		blackBackground = geDi.getNextBoolean();

		double actPercentValue = thresholdValuePercentual;
		double actAbsoluteValue = thresholdValueAbsolute;
	
		Choice modusChoice = (Choice) geDi.getChoices().get(0);
		Scrollbar absScroll = (Scrollbar) geDi.getSliders().get(1);
		TextField absTextField = (TextField) geDi.getNumericFields().get(1);
		Scrollbar relScroll = (Scrollbar) geDi.getSliders().get(0);
		TextField relTextField = (TextField) geDi.getNumericFields().get(0);

		boolean relSelected = (modusChoice.getSelectedItem()==modusChoice.getItem(0));
		absScroll.setEnabled(!relSelected);
		absTextField.setEnabled(!relSelected);
		relScroll.setEnabled(relSelected);
		relTextField.setEnabled(relSelected);
		

	
		//(Choice)(geDi.getChoices().get(0))
		// Validierung
		if (Double.isNaN(actPercentValue) || actPercentValue < 0) {
			IJ.showMessage("'Relativ proportion have to be > 0!");
			return false;
		}
		
		if (Double.isNaN(actAbsoluteValue) || actAbsoluteValue < 0) {
			IJ.showMessage("'Absolute number have to be > 0!");
			return false;
		}
		
		return true;
	
	}

	@Override
	public void setNPasses(int nPasses) {
		// TODO Auto-generated method stub
		
	}
	
	
	private void doFourierFilter(ImageProcessor ip) {
		
		shapeSmoothingUtil.setDrawOnlyContours(drawOnlyContours);
		shapeSmoothingUtil.setBlackBackground(blackBackground);
		
		if (doAbsoluteThreshold) {
			shapeSmoothingUtil.fourierFilter(ip, thresholdValueAbsolute, false, doOutputDescriptors && !previewing);
		} else {
			shapeSmoothingUtil.fourierFilter(ip, thresholdValuePercentual, true, doOutputDescriptors && !previewing);
		}

	}
	
	/**
	 * Dupliziert ein Bild (die Kopie erscheint in einem neuen Fenster)
	 * 
	 * @param sourceIp {@link ImageProcessor} des zu duplizierenden Bildes
	 * @param newTitle Überschrift für das neue Bildfenster
	 * @return {@link ImagePlus} des Duplikats
	 */
	public ImagePlus duplicateWindow (ImageProcessor sourceIp, String newTitle) {
		ImagePlus newImp = NewImage.createByteImage(newTitle, sourceIp.getWidth(), sourceIp.getHeight(), 1, NewImage.FILL_BLACK);
		newImp.setProcessor(sourceIp.duplicate());		
		return newImp;
	}

	

	
}