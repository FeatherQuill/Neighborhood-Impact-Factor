
public class Cell {
	/** This should be the stain which we think influences the neighborhood (GATA6, HA, etc). */
	private double independentMarkerIntensity;
	/** This should be the stain for one of the cell fates (CEBPa, etc). */
	private double marker1Intensity;
	/** This should be the stain for the other cell fate (CD34, etc). */
	private double marker2Intensity;
	private double x;
	private double y;
	private int objNumber;
	private int imageNumber;
	
	private double if1 = -1;
	private double if2 = -1;
	private double if3 = -1;
	private double if4 = -1;
	
	public Cell(int id, int imageId, double x, double y, double independentVal, double marker1Val, double marker2Val){
		objNumber = id;
		imageNumber = imageId;
		this.x = x;
		this.y = y;
		independentMarkerIntensity = independentVal;
		marker1Intensity = marker1Val;
		marker2Intensity = marker2Val;
	}
	
	/**
	 * @return The intensity of the expression of the control gene
	 * strain (such as GATA6, or HA).
	 */
	public double getIndependentStainExpression(){
		return independentMarkerIntensity;
	}
	
	/**
	 * @return The intensity of the stain for one of the cells fates
	 * (CEBPa, etc).
	 */
	public double getMarker1Intensity(){
		return marker1Intensity;
	}
	
	/**
	 * @return The intensity of the stain for the other cells fates
	 * (CD34, etc).
	 */
	public double getMarker2Intensity(){
		return marker2Intensity;
	}
	
	public double getX(){
		return x;
	}
	
	public double getY(){
		return y;
	}
	
	public int getID(){
		return objNumber;
	}
	
	public String toString(){
		return "ObjNum: "+objNumber+"\nX,Y: "+x+", "+y+"\nIndependent Marker: "+independentMarkerIntensity+"\nMarker 1: "+marker1Intensity+"\nMarker 2: "+marker2Intensity;
	}
	
	public boolean isMarker1Pos(){
		return marker1Intensity>0;
	}
	
	public boolean isMarker2Pos(){
		return marker2Intensity>0;
	}
	
	public boolean isIndCellPos(){
		return independentMarkerIntensity>0;
	}
	
	public void setIFValue(double val, String ifType){
		if(ifType.equals(Main.IF1)){
			if1 = val;
			return;
		}
		if(ifType.equals(Main.IF2)){
			if2 = val;
			return;
		}
		if(ifType.equals(Main.IF3)){
			if3 = val;
			return;
		}
		if(ifType.equals(Main.IF4)){
			if4 = val;
			return;
		}
	}
	
	public double getIF(String ifType){
		if(ifType.equals(Main.IF1)){
			return if1;
		}
		if(ifType.equals(Main.IF2)){
			return if2;
		}
		if(ifType.equals(Main.IF3)){
			return if3;
		}
		if(ifType.equals(Main.IF4)){
			return if4;
		}
		return -1;
	}
	
	public int getImageId(){
		return imageNumber;
	}
}
