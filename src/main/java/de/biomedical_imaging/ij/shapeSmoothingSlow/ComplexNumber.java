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
