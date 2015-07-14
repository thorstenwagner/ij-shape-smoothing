package de.biomedical_imaging.ij.shapeSmoothing.tests;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import ij.IJ;
import ij.ImagePlus;
import ij.blob.Blob;
import ij.blob.ManyBlobs;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import org.junit.Before;
import org.junit.Test;

import de.biomedical_imaging.ij.shapeSmoothing.ShapeSmoothingUtil;
import de.biomedical_imaging.ij.shapeSmoothingPlugin.Shape_Smoothing;

public class ShapeSmoothingUtilTest {
	
	private ImagePlus testImp;
	private ImageProcessor testIp;
	private int width;
	private int height;
	private ImageProcessor newIp;
	private ShapeSmoothingUtil shapeSmoothingUtil;
	private Shape_Smoothing shapeSmoothing;
	
	public ShapeSmoothingUtilTest() {
		testImp = IJ.openImage(ShapeSmoothingUtilTest.class.getResource("testPic.tif").getPath());
		testIp = testImp.getProcessor();
		width = testIp.getWidth();
		height = testIp.getHeight();
		shapeSmoothingUtil = new ShapeSmoothingUtil();
		shapeSmoothing = new Shape_Smoothing();
	}
	
	/**
	 * Prüft, ob zwei Bilder (bzw. zwei Pixel-Arrays) gleich sind
	 * 
	 * @param image1 erstes Bild
	 * @param image2 zweites Bild
	 * @return true, wenn die Abmessungen der beiden Bilder gleich sind und alle Pixeln gleiche Werte haben, sonst false
	 */
	private boolean testImageEquality(int[][] image1, int[][] image2) {
		// Abmessungen prüfen
		if (image1.length != image2.length || image1[0].length != image2[0].length) {
			return false;
		}
		
		// Pixelwerte prüfen
		int unequalPxls = 0;
		for (int x = 0; x < image1.length; x++) {
			for (int y = 0; y < image1[0].length; y++) {
				if (image1[x][y] != image2[x][y]) {
					System.out.println("Ungleicher Pixel: " + x + ", " + y);
					unequalPxls++; //return false;
				}
			}
		}
		if (unequalPxls > 0) {
			System.out.println("Anzahl ungleicher Pixel: " + unequalPxls);
			return false;
		}
		
		return true;
	}

	@Before
	public void setUp() throws Exception {		
		newIp = new ByteProcessor(testImp.getWidth(), testImp.getHeight());
		newIp.invert(); // invertieren, da der Konstruktor ein schwarzes Bild liefert
	}
	
	@Test
	public void fourierFilter_PointTest() {		
		
		
		// Prüfe, ob die Pixeln auf den Koordinaten der Schwerpunkte der Formen schwarz sind
		ManyBlobs allBlobs = new ManyBlobs(testImp);
		allBlobs.findConnectedComponents();
		int sumXCoord = 0;
		int sumYCoord = 0;
		for(Blob blob: allBlobs) {
			int centerX = (int) Math.round(blob.getCenterOfGravity().getX());
			int centerY = (int) Math.round(blob.getCenterOfGravity().getY());
			sumXCoord +=centerX;
			sumYCoord += centerY;
		}
		//shapeSmoothingUtil.fourierFilter(testImp, newIp, 1, false);
		shapeSmoothingUtil.fourierFilter(testImp.getProcessor(), 1, false, false);
		int sumFourierXCoord = 0;
		int sumFourierYCoord = 0;
		for (int x = 0; x < width; x++) {
			for (int y = 0; y < height; y++) {
				if (testImp.getProcessor().getPixel(x, y) == 0) {
					sumFourierXCoord += x;
					sumFourierYCoord += y;
				}
			}
		}
		int diff = Math.abs(sumFourierXCoord-sumXCoord) + Math.abs(sumFourierYCoord-sumYCoord);
		//Eine Differenz entsteht dadurch, dass die Konturen vor der Fouier Transformation noch transformiert
		//werden (Äquidistanz). Zur Berechnung des Schwerpunkts über IJBlob aber nicht.
		assertEquals(0, diff, 6);
	}
	
	@Test
	public void fourierFilter_CircleTest() {
		ManyBlobs allBlobsOriginal = new ManyBlobs(testImp);
		allBlobsOriginal.findConnectedComponents();
		
	
		shapeSmoothingUtil.fourierFilter(testImp.getProcessor(), 2, false, false);
		
		ManyBlobs allBlobsCircles = new ManyBlobs(testImp);
		allBlobsCircles.findConnectedComponents();
		
		// Prüfen, dass es genauso viele Kreise gibt, wie Figuren
		assertEquals(allBlobsOriginal.size(), allBlobsCircles.size());
		
		for (int i = 0; i < allBlobsCircles.size(); i++) {
			Blob circle = allBlobsCircles.get(i);
			assertEquals(1, circle.getThinnesRatio(), 0.1); // Prüfen, dass es sich tatsächlich um Kreise handelt
			circle.getCenterOfGravity().equals(allBlobsOriginal.get(i).getCenterOfGravity()); // Prüfen, ob Schwerpunkte jeweils gleich sind
		}
	}
	
	@Test
	public void fourierFilter_AllFDsTest() {
		ManyBlobs allBlobs = new ManyBlobs(testImp);
		allBlobs.findConnectedComponents();	
		
		
		
		
		//Das Ergebnis sollte nur mit den Konturen der Originalformen verglichen werden, um Unstimmigkeiten zu vermeiden
		////shapeSmoothing.duplicateWindow(testIp, "Konturen");
		ImageProcessor newIp = new ByteProcessor(testImp.getWidth(), testImp.getHeight());	
		newIp.invert();
		
		for (Blob blob: allBlobs) {
			newIp.drawPolygon(shapeSmoothingUtil.toEquidistantPolygon(blob.getOuterContour()));
		}
		ImagePlus newImp = new ImagePlus("", newIp);
		
		// Die Konturen an sich sind nicht (immer) skeletiert
		//IJ.run(newImp, "Skeletonize", "");
	
		int maxNumOfContourPoints = 0;
		for (Blob blob: allBlobs) {
			int blobsNumOfContourPoints = blob.getOuterContour().npoints;
			if (blobsNumOfContourPoints > maxNumOfContourPoints) {
				maxNumOfContourPoints = blobsNumOfContourPoints;
			}			
		}
		shapeSmoothingUtil.setDrawOnlyContours(true);
		shapeSmoothingUtil.fourierFilter(testImp.getProcessor(), (double)maxNumOfContourPoints, false, false);

		assertTrue(testImageEquality(testImp.getProcessor().getIntArray(), newIp.getIntArray()));
	}

}
