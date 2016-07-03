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

package de.biomedical_imaging.ij.shapeSmoothingSlow;
import ij.measure.ResultsTable;
import ij.plugin.filter.Analyzer;

import java.util.Vector;


public class MyUsefulMethods {
	
	// DFT-Hintransformation
	public static Vector<ComplexNumber> dft (Vector<ComplexNumber> contourPoints) {
		if (contourPoints == null || contourPoints.size() == 0) {
			return null;
		}
		
		Vector<ComplexNumber> fds = new Vector<ComplexNumber>();
		int N = contourPoints.size();
		
		ComplexNumber fd = null;
		ComplexNumber summand = null;
		double angle = 0;
		for (int k = 0; k < N; k++) {
			fd = new ComplexNumber();
			for (int j = 0; j < N; j++) {
				summand = contourPoints.get(j);
				angle = (2.0*Math.PI*k*j) / (1.0 *N);
				summand = summand.multiplication(
							new ComplexNumber(Math.cos(angle), -1.0 * Math.sin(angle))
						);
				fd = fd.addition(summand);
			}
			//fd = fd.multiplication(1.0 / N);
			fds.add(fd);
		}
		return fds;
	}
	
	// DFT-Ruecktransformation
	public static Vector<ComplexNumber> dftInverse (Vector<ComplexNumber> fds) {
		if (fds == null || fds.size() == 0) {
			return null;
		}
		
		Vector<ComplexNumber> contourPoints = new Vector<ComplexNumber>();
		int N = fds.size();
		
		ComplexNumber point = null;
		ComplexNumber summand = null;
		double angle = 0;
		for (int k = 0; k < N; k++) {
			point = new ComplexNumber();
			for (int j = 0; j < N; j++) {
				summand = fds.get(j);
				angle = (2.0*Math.PI*k*j) / (1.0 *N);
				summand = summand.multiplication(
							new ComplexNumber(Math.cos(angle), Math.sin(angle))
						);
				point = point.addition(summand);
			}
			point = point.multiplication(1.0/N);
			contourPoints.add(point);
		}
		return contourPoints;
	}
	
	public static void showComplexVectorInTable (Vector<ComplexNumber> vector, String tableTitle) {
		ResultsTable rt = new ResultsTable();
		Analyzer.setResultsTable(rt);
		
		for (ComplexNumber cn: vector) {
			rt.incrementCounter();
			rt.addValue("Reell (x)", cn.getRe());
			rt.addValue("Imaginär (y)", cn.getIm());			
		}
		rt.show(tableTitle);
	}
		
}
