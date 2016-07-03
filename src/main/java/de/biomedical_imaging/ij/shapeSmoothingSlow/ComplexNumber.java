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
public class ComplexNumber{
	
	private double re;		
	private double im;
	
	private static final double SIGMA = 0.00001;
	
	public ComplexNumber(double re, double im) {
		this.re = re;
		this.im = im;
	}
	
	public ComplexNumber(){
		this(0.0, 0.0);
	}
	
	public double getRe() {
		return re;
	}		
	public double getIm() {
		return im;
	}
	
	public ComplexNumber addition (ComplexNumber summand) {
		ComplexNumber result = new ComplexNumber();
		result.re = this.re + summand.re;
		result.im = this.im + summand.im;			
		return result;
	}
	
	public ComplexNumber multiplication (ComplexNumber factor) {
		ComplexNumber result = new ComplexNumber();
		result.re = this.re*factor.re - this.im*factor.im;
		result.im = this.re*factor.im + this.im*factor.re;			
		return result;
	}
	
	public ComplexNumber multiplication (double factor) {
		ComplexNumber result = new ComplexNumber();
		result.re = factor*this.re;
		result.im = factor*this.im;			
		return result;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj instanceof ComplexNumber) {
			ComplexNumber other = (ComplexNumber) obj;
			return (Math.abs(other.re-this.re) < SIGMA && Math.abs(other.im-this.im) < SIGMA);
		} else {
			return false;
		}
	}
	
	public String toString() {
		String plusSign = new String();
		if (im >= 0) {
			plusSign = "+";
		}
		return re + plusSign + im + "i";
	}
	
	public double getAbsVal() {
		return Math.sqrt(Math.pow(re, 2) + Math.pow(im, 2));		
	}
	
	public boolean greaterThen(ComplexNumber otherCN){
		return this.getAbsVal() > otherCN.getAbsVal();			
	}
	
}
