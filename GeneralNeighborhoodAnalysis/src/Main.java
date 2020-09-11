import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.text.DecimalFormat;
import java.util.ArrayList;

import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RVector;
import org.rosuda.JRI.Rengine;
import org.rosuda.REngine.REngine;

public class Main {
	private static String dbPath;
	private static String dbURL;
	private static Connection c;
	private static String query;
	
	private static String if1Output;
	private static String if2Output;
	private static String if3Output;
	private static String if4Output;
	private static String csvHeader;
	
	protected static final String IF1 = "Total Expression";
	protected static final String IF2 = "Local Density";
	protected static final String IF3 = "Distance Reversed Expression";
	protected static final String IF4 = "Distance Adjusted Expression";
	
	private static String MARKER1; //= "cebpa+";
	private static String MARKER2; //= "cd34+";
//	private static final String MEAN_HA_BOTH = "both";
	
	private static String if1BoundaryOutput;
	private static String if2BoundaryOutput;
	private static String if3BoundaryOutput;
	private static String if4BoundaryOutput;
	private static String boundarySummaryHeader;
	private static String boundarySummaryOutput;
	
	private static ArrayList<Double> radiiForBoundary; 
	private static String boundaryHeader;
	
	private static ArrayList<Cell> allCells;
	private static int indCells;
	private static ArrayList<Cell> popIndCells;
	
	private static ArrayList<Cell> popIndCellsMarker1;
	private static ArrayList<Cell> popIndCellsMarker2;
	private static ArrayList<Cell> popIndCellsBoth;
	private static ArrayList<Cell> popIndCellsNone;
	
	private static String summaryOutput;
	
	private static double indMean;
	private static double indMax;
	private static double indMeanLow;
	private static double indMeanHigh;
	private static double populationWidth;
	private static ArrayList<Double> indMeans;
	private static ArrayList<Double> indMeanLows;
	private static ArrayList<Double> indMeanHighs;
	private static String populationHeader;
	private static String populationOutput;
	
	private static String plotCsvOutput;
	
	private static String indCellType; // = "nuclei_children_ha_stain_count";
	private static String indType;
	private static String indIntensity; // = "nuclei_math_normalized_ha";
	
	//TODO: Make sure you set this before you run the analysis
	private static double neighborhoodRadius;// = 770.0;
	
	private static Rengine rengine;
	
	private static double m1If1Avg = 0;
	private static double m2If1Avg = 0;
	private static double m1If2Avg = 0;
	private static double m2If2Avg = 0;
	private static double m1If3Avg = 0;
	private static double m2If3Avg = 0;
	private static double m1If4Avg = 0;
	private static double m2If4Avg = 0;
	
	
	private static int neighborCount = 0;
	
	public static void main(String[] args){
		new MainFrame();
	}
	
	
	
	private void oldMain(){
		//TODO: make sure these are set correctly
				String basePath = "C:/Users/shay_/Documents/ASU/thesis/Images for Analysis/thesis stainings/11 - PGP1 D3 - rbT, msFOXA2, gtGATA6/cropped/standard_crop2/output/";
//				String basePath = "C:/Users/shay_/Documents/ASU/thesis/Images for Analysis/Early Time Points/Day 3/PGP1 D3, T, GATA6, FOXA2, reimaged 11-11-16/cropped/standard_crop2/output/";
//				String basePath = "C:/Users/shay_/Documents/ASU/thesis/Images for Analysis/CEBPa-CD34-HA-DAPI/cropped/standard_crop2/output/";
//				dbPath = basePath+"pgp1-day3-t-foxa2-gata6-standard_crop2.db";
//				dbPath = basePath+"standard_crop2.db";
				
				
				String outputPath = basePath+"NIF/";
				//make sure the output directory exists
				new File(outputPath).mkdirs();
				
//				dbPath = basePath+"standard_crop.db";
				dbURL = "jdbc:sqlite:"+dbPath;
				
				plotCsvOutput = outputPath+"plot.csv";
				
				
				//TODO: set these before the run
				MARKER1 = "T";
				MARKER2 = "FOXA2";
				indType = "GATA6";
				indCellType = "nuclei_children_"+indType+"_stain_count";
				indIntensity = "nuclei_math_normalized_"+indType;
				
				neighborhoodRadius = 30;
				
				//----------------------------------------------//
				
				
				radiiForBoundary = new ArrayList<Double>();
				//populate a radius list for boundary testing
				for(int i=1; i<282; i++){
					double number = i*10;
					radiiForBoundary.add(number);
				}
//				radiiForBoundary.add(10.0);
//				radiiForBoundary.add(20.0);
//				radiiForBoundary.add(30.0);
//				radiiForBoundary.add(40.0);
//				radiiForBoundary.add(50.0);
				
				

				
				boundarySummaryHeader = "IF Method, Max % Diff, Radius";
				boundarySummaryOutput = outputPath+"Boundary_Summary.txt";
				
				popIndCellsBoth = new ArrayList<Cell>();
				popIndCellsMarker1 = new ArrayList<Cell>();
				popIndCellsMarker2 = new ArrayList<Cell>();
				popIndCellsNone = new ArrayList<Cell>();
				
				allCells = getAllCells();
				indCells = getNumberOfCellType(indCellType, indIntensity);
				indMean = getMeanInd();
				//TODO: if you don't want to use the mean HA, override it here
//				haMean = 0.183;
//				indMean = 0.15;
				populationWidth = indMax*.005;
				System.out.println("Pop width: "+populationWidth);
				indMeanLow = indMean-populationWidth/2;
				indMeanHigh = indMean+populationWidth/2;
				popIndCells = getPopIndCells(indMeanLow, indMeanHigh);
				
				indMeanLows = new ArrayList<Double>();
				indMeanHighs = new ArrayList<Double>();
				indMeans = new ArrayList<Double>();
				//step through min to max ha in 5% increments
				System.out.println("Max "+indType+": "+indMax);
				double step = indMax*.05;
				System.out.println("Step: "+step);
				for(double i=0; i<indMax; i+=step){
					double mean = i;
//					System.out.println("mean: "+mean);
					double low;
					double high;
					//set low
					if(i>step/2){
						low = mean-step/2;
					}else{
						low = 0;
					}
					//set high
					if(i<mean+step/2){
						high = mean+step/2;
					}else{
						high = indMax;
					}
					
					//add them to lists
					indMeanLows.add(low);
					indMeans.add(mean);
					indMeanHighs.add(high);
				}
				populationHeader = "Center "+indType+", "+MARKER1+"+ Cells, "+MARKER2+"+ Cells";
				populationOutput = outputPath+"Population_Sample_Comparison.csv";

				DecimalFormat df = new DecimalFormat("#0.0000");

				
				boundaryHeader = "Radius (pixels), Ave "+MARKER1+"+ IF, Ave "+MARKER2+"+ IF, % Difference, P-Value, Ave # of Cells in Neighborhood";
				if1BoundaryOutput = outputPath+"Total_Expression_IF_Boundary_"+indType+"-"+df.format(indMean)+".csv";
				if2BoundaryOutput = outputPath+"Local_Density_IF_Boundary_"+indType+"-"+df.format(indMean)+".csv";
				if3BoundaryOutput = outputPath+"Distance_Reversed_Expression_IF_Boundary_"+indType+"-"+df.format(indMean)+".csv";
				if4BoundaryOutput = outputPath+"Distance_Adjusted_Expression_IF_Boundary_"+indType+"-"+df.format(indMean)+".csv";
				
				
				if1Output = outputPath+"Total_Expression_IF_"+indType+"-"+df.format(indMean)+"_rad-"+neighborhoodRadius+".csv";
				if2Output = outputPath+"Local_Density_IF_"+indType+"-"+df.format(indMean)+"_rad-"+neighborhoodRadius+".csv";
				if3Output = outputPath+"Distance_Reversed_Expression_IF_"+indType+"-"+df.format(indMean)+"_rad-"+neighborhoodRadius+".csv";
				if4Output = outputPath+"Distance_Adjusted_Expression_IF_"+indType+"-"+df.format(indMean)+"_rad-"+neighborhoodRadius+".csv";
				csvHeader = MARKER1+"+ "+indType+", IF Value, "
							+MARKER2+"+ "+indType+", IF Value, "
							+MARKER1+"+/"+MARKER2+"+ "+indType+", IF Value, "
							+MARKER1+"-/"+MARKER2+"- "+indType+", IF Value";
				
				System.out.println("Total Cells: "+allCells.size());
				System.out.println(indType+"+ Cells: "+indCells);
				System.out.println("Mean "+indType+": "+indMean);

				System.out.println("Mean "+indType+" Cells: "+popIndCells.size());
				System.out.println("Mean "+indType+" Pure "+MARKER1+"+ Cells: "+popIndCellsMarker1.size());
				System.out.println("Mean "+indType+" Pure "+MARKER2+"+ Cells: "+popIndCellsMarker2.size());
				System.out.println("Mean "+indType+" Mixed "+MARKER1+"+/"+MARKER2+"+ Cells: "+popIndCellsBoth.size());
				System.out.println("Mean "+indType+" Other Cells: "+popIndCellsNone.size());
				
				summaryOutput = outputPath+"Population_Summary.txt";
				
				try {
					writePlotCsv();
					writeIFFile(IF1);
					writeIFFile(IF2);
					writeIFFile(IF3);
					writeIFFile(IF4);
					writeSummaryFile();
//					writeIFBoundaryFile(IF1);
//					writeIFBoundaryFile(IF2);
//					writeIFBoundaryFile(IF3);
//					writeIFBoundaryFile(IF4);
//					writePopulationComparisonFile();
//					writeBoundarySummaryFile();
//					rengine.end();
				} catch (IOException e) {
					e.printStackTrace();
				}
				
				
//				System.out.println("Total Expression");
//				System.out.println(MARKER1+": "+getAveIFForPop(popIndCellsMarker1, IF1));
//				System.out.println(MARKER2+": "+getAveIFForPop(popIndCellsMarker2, IF1));
	}
	
	
	private static void writePlotCsv() throws IOException{
		FileWriter fw1 = new FileWriter(plotCsvOutput);
		BufferedWriter bw1 = new BufferedWriter(fw1);
		
		//write header
		bw1.write(MARKER1+","+MARKER2+","+indType+","+MARKER1+","+indType+","+MARKER2);
		bw1.newLine();
		for(Cell c : allCells){
			//marker1 and marker2
			if(c.isMarker1Pos() && c.isMarker2Pos()){
				bw1.write(c.getMarker1Intensity()+","+c.getMarker2Intensity()+",");
			}else{
				bw1.write(",,");
			}
			//marker1 and indCell
			if(c.isMarker1Pos() && c.isIndCellPos()){
				bw1.write(c.getMarker1Intensity()+","+c.getIndependentStainExpression()+",");
			}else{
				bw1.write(",,");
			}
			//marker2 and indcell
			if(c.isMarker2Pos() && c.isIndCellPos()){
				bw1.write(c.getMarker2Intensity()+","+c.getIndependentStainExpression());
			}else{
				bw1.write(",");
			}
			bw1.newLine();
		}
		
		bw1.close();
		fw1.close();
		System.out.println("Plot file written successfully!");
	}
	
	
	private static double getAveIFForPop(ArrayList<Cell> pop, String ifType){
		double ave = 0;
		double total = 0;
		
		for(Cell c : pop){
			total += c.getIF(ifType);
		}
		
		ave = total/pop.size();
		
		return ave;
	}
	
	
	private static void writeSummaryFile() throws IOException{
		
		//Calculate the average number of neighbors per marker1+ and marker2+ cells
		int m1Size = popIndCellsMarker1.size();
		int m2Size = popIndCellsMarker2.size();
		int max = Math.max(m1Size, m2Size);
		//set count to 0 in case these methods have already been run
		neighborCount = 0;
		for (int i=0; i<max; i++){
			if(m1Size>i){
				setIFOfTargetCell(IF1, popIndCellsMarker1.get(i), neighborhoodRadius);
			}
			if(m2Size>i){
				setIFOfTargetCell(IF1, popIndCellsMarker2.get(i), neighborhoodRadius);
			}
		}
		
		//find average
		double avgNeighbors = ((double)neighborCount)/(m1Size+m2Size);
		
		//find p values for each IF method
		ArrayList pVals = calculatePValues();
		
		FileWriter fw = new FileWriter(summaryOutput);
		BufferedWriter bw = new BufferedWriter(fw);
		
		bw.write("Total Cells: "+allCells.size());
		bw.newLine();
		bw.write(indType+"+ Cells: "+indCells);
		bw.newLine();
		bw.write("Mean "+indType+": "+indMean);
		bw.newLine();
		bw.write("GATA6 max: "+indMax);
		bw.newLine();
		bw.write("Pop width: "+populationWidth);
		bw.newLine();
		bw.write("Mean "+indType+" Cells: "+popIndCells.size());
		bw.newLine();
		bw.write("Mean "+indType+" Pure "+MARKER1+"+ Cells: "+popIndCellsMarker1.size());
		bw.newLine();
		bw.write("Mean "+indType+" Pure "+MARKER2+"+ Cells: "+popIndCellsMarker2.size());
		bw.newLine();
		bw.write("Mean "+indType+" Mixed "+MARKER1+"+/"+MARKER2+"+ Cells: "+popIndCellsBoth.size());
		bw.newLine();
		bw.write("Mean "+indType+" Other Cells: "+popIndCellsNone.size());
		bw.newLine();
		bw.write("Neighborhood size (pixels): "+neighborhoodRadius); 
		bw.newLine();
		bw.write("Average number of neightbors: "+avgNeighbors);
		bw.newLine();
		bw.newLine();
		bw.write(MARKER1+" mean IF1: "+m1If1Avg);
		bw.newLine();
		bw.write(MARKER2+" mean IF1: "+m2If1Avg);
		bw.newLine();
		bw.write("IF1 p-value: "+pVals.get(0));
		bw.newLine();
		bw.newLine();
		bw.write(MARKER1+" mean IF2: "+m1If2Avg);
		bw.newLine();
		bw.write(MARKER2+" mean IF2: "+m2If2Avg);
		bw.newLine();
		bw.write("IF2 p-value: "+pVals.get(1));
		bw.newLine();
		bw.newLine();
		bw.write(MARKER1+" mean IF3: "+m1If3Avg);
		bw.newLine();
		bw.write(MARKER2+" mean IF3: "+m2If3Avg);
		bw.newLine();
		bw.write("IF3 p-value: "+pVals.get(2));
		bw.newLine();
		bw.newLine();
		bw.write(MARKER1+" mean IF4: "+m1If4Avg);
		bw.newLine();
		bw.write(MARKER2+" mean IF4: "+m2If4Avg);
		bw.newLine();
		bw.write("IF4 p-value: "+pVals.get(3));
		bw.newLine();
		bw.newLine();
		
		bw.close();
		fw.close();
		System.out.println("Summary File Complete!");
	}
	
	
	public static ArrayList<Double> calculatePValues(){
		ArrayList<Double> result = new ArrayList<Double>();
		
		//create double[] for each IF value for each cell type
		double[] marker1IF1 = new double[popIndCellsMarker1.size()];
		double[] marker1IF2 = new double[popIndCellsMarker1.size()];
		double[] marker1IF3 = new double[popIndCellsMarker1.size()];
		double[] marker1IF4 = new double[popIndCellsMarker1.size()];
		double[] marker2IF1 = new double[popIndCellsMarker2.size()];
		double[] marker2IF2 = new double[popIndCellsMarker2.size()];
		double[] marker2IF3 = new double[popIndCellsMarker2.size()];
		double[] marker2IF4 = new double[popIndCellsMarker2.size()];
		
		m1If1Avg = 0;
		m2If1Avg = 0;
		m1If2Avg = 0;
		m2If2Avg = 0;
		m1If3Avg = 0;
		m2If3Avg = 0;
		m1If4Avg = 0;
		m2If4Avg = 0;
		//populate arrays
		int i=0;
		for(Cell c : popIndCellsMarker1){
			marker1IF1[i] = c.getIF(IF1);
			marker1IF2[i] = c.getIF(IF2);
			marker1IF3[i] = c.getIF(IF3);
			marker1IF4[i] = c.getIF(IF4);
			
			//add totals so avg can be calculated
			m1If1Avg += c.getIF(IF1);
			m1If2Avg += c.getIF(IF2);
			m1If3Avg += c.getIF(IF3);
			m1If4Avg += c.getIF(IF4);
			
			i++;
		}
		i=0;
		for(Cell c: popIndCellsMarker2){
			marker2IF1[i] = c.getIF(IF1);
			marker2IF2[i] = c.getIF(IF2);
			marker2IF3[i] = c.getIF(IF3);
			marker2IF4[i] = c.getIF(IF4);
			
			//add totals so avg can be calculated
			m2If1Avg += c.getIF(IF1);
			m2If2Avg += c.getIF(IF2);
			m2If3Avg += c.getIF(IF3);
			m2If4Avg += c.getIF(IF4);
			
			i++;
		}		
		
		m1If1Avg = m1If1Avg/popIndCellsMarker1.size();
		m1If2Avg = m1If2Avg/popIndCellsMarker1.size();
		m1If3Avg = m1If3Avg/popIndCellsMarker1.size();
		m1If4Avg = m1If4Avg/popIndCellsMarker1.size();
		m2If1Avg = m2If1Avg/popIndCellsMarker2.size();
		m2If2Avg = m2If2Avg/popIndCellsMarker2.size();
		m2If3Avg = m2If3Avg/popIndCellsMarker2.size();
		m2If4Avg = m2If4Avg/popIndCellsMarker2.size();
		
		System.out.println("setting up r engine");
		Rengine engine = new Rengine(new String[] {"--no-save"}, false, null);

		//repeat the assignment and evaluation 4 times (one for each IF method)
		for(int j=0; j<4; j++){
			//IF1
			if(j==0){
				engine.assign("a1", marker1IF1);
				engine.assign("a2", marker2IF1);
				System.out.println("m1If1Avg: "+m1If1Avg);
				System.out.println("m2If1Avg: "+m2If1Avg);
			}
			//IF2
			if(j==1){
				engine.assign("a1", marker1IF2);
				engine.assign("a2", marker2IF2);
				System.out.println("m1If2Avg: "+m1If2Avg);
				System.out.println("m2If2Avg: "+m2If2Avg);
			}
			//IF3
			if(j==2){
				engine.assign("a1", marker1IF3);
				engine.assign("a2", marker2IF3);
				System.out.println("m1If3Avg: "+m1If3Avg);
				System.out.println("m2If3Avg: "+m2If3Avg);
			}
			//IF4
			if(j==3){
				engine.assign("a1", marker1IF4);
				engine.assign("a2", marker2IF4);
				System.out.println("m1If4Avg: "+m1If4Avg);
				System.out.println("m2If4Avg: "+m2If4Avg);
			}
			
			REXP rexp = engine.eval("t.test(a1, a2, paired=FALSE)");
			RVector rvec = rexp.asVector();
			
			//p value is the 3rd value
			String pString = rvec.get(2).toString();
			//Take off beginning "[REAL* (" and end: ")]"
			pString = pString.substring(8, pString.length()-2);
			System.out.println("pval: "+pString);
			result.add(Double.parseDouble(pString));
		}
		
		engine.end();
		
		return result;
	}
	
	
	private static double calculatePValue(String ifType){
		double result = 0;
		
		double[] marker1IF = new double[popIndCellsMarker1.size()];
		double[] marker2IF = new double[popIndCellsMarker2.size()];
		
		System.out.println("size of marker1: "+marker1IF.length);
		System.out.println("size of marker2: "+marker2IF.length);
		
		//populate arrays
		int i=0;
		for(Cell c : popIndCellsMarker1){
			marker1IF[i] = c.getIF(ifType);
			i++;
		}
		i=0;
		for(Cell c : popIndCellsMarker2){
			marker2IF[i] = c.getIF(ifType);
			i++;
		}
		
		if(rengine==null){
			System.out.println("setting up r engine");
			rengine = new Rengine(new String[] {"--no-save"}, false, null);
		}
		
		rengine.assign("a1", marker1IF);
		rengine.assign("a2", marker2IF);
		
		REXP rexp = rengine.eval("t.test(a1, a2, paired=FALSE)");
		RVector rvec = rexp.asVector();
		
//		System.out.println("rexp: "+rvec.toString());
		
		//p value is the 3rd value
		String pString = rvec.get(2).toString();
		//Take off beginning "[REAL* (" and end: ")]"
		pString = pString.substring(8, pString.length()-2);
		System.out.println("pval: "+pString);
		result = Double.parseDouble(pString);
		
		
		
		return result;
	}
	
	
		
	/**
	 * Writes out a csv with the following header:
	 * "CEBPa+ HA, IF Value, CD34+ HA, IF Value, CEBPa+/CD34+ HA,
	 * IF Value, CEBPa-/CD34- HA, IF Value".
	 * IF1 is calculated by using the Total Expression 
	 * method.
	 * @throws IOException
	 */
	private static void writeIFFile(String ifType) throws IOException{
		String fileName;
		switch(ifType){
			case IF1: 	fileName = if1Output;
						break;
			case IF2: 	fileName = if2Output;
						break;
			case IF3:	fileName = if3Output;
						break;
			case IF4:	fileName = if4Output;
						break;
			default:	fileName = "";
		}
		FileWriter fw = new FileWriter(fileName);
		BufferedWriter bw = new BufferedWriter(fw);
		//write header
		bw.write(csvHeader);
		bw.newLine();
		int maxIndex = Math.max(Math.max(popIndCellsBoth.size(), popIndCellsNone.size()), Math.max(popIndCellsMarker1.size(),	popIndCellsMarker2.size()));
		
		for(int i=0; i<maxIndex; i++){
			//if there are marker 1+ cells left, add to csv
			if(popIndCellsMarker1.size()>i){
				//get cell
				Cell cell = popIndCellsMarker1.get(i);
				//find the IF for the cell
				setIFOfTargetCell(ifType, cell, neighborhoodRadius);
				//print out the cell expression and it's IF
				bw.write(cell.getIndependentStainExpression()+","+cell.getIF(ifType));
			}else{
				bw.write(",");
			}
			//add marker 2+ cells if there are any
			if(popIndCellsMarker2.size()>i){
				//get cell
				Cell cell = popIndCellsMarker2.get(i);
				//set IF for the cell
				setIFOfTargetCell(ifType, cell, neighborhoodRadius);
				//print expression and IF
				bw.write(","+cell.getIndependentStainExpression()+","+cell.getIF(ifType));
			}else{
				bw.write(", ,");
			}
			//add both positive cells if there are any
			if(popIndCellsBoth.size()>i){
				//get cell
				Cell cell = popIndCellsBoth.get(i);
				//set IF for the cell
				setIFOfTargetCell(ifType, cell, neighborhoodRadius);
				//print expression and IF
				bw.write(","+cell.getIndependentStainExpression()+","+cell.getIF(ifType));
			}else{
				bw.write(", ,");
			}
			//add both negative cells if there are any
			if(popIndCellsNone.size()>i){
				//get cell
				Cell cell = popIndCellsNone.get(i);
				//set IF for cell
				setIFOfTargetCell(ifType, cell, neighborhoodRadius);
				//print expression and IF
				bw.write(","+cell.getIndependentStainExpression()+","+cell.getIF(ifType));
			}else{
				bw.write(", ,");
			}
			bw.newLine();
		}
		
		bw.close();
		fw.close();
		System.out.println(fileName+" complete!");
	}
	
	
	private static void writeIFBoundaryFile(String ifType) throws IOException{
		String fileName;
		switch(ifType){
			case IF1: 	fileName = if1BoundaryOutput;
						break;
			case IF2: 	fileName = if2BoundaryOutput;
						break;
			case IF3:	fileName = if3BoundaryOutput;
						break;
			case IF4:	fileName = if4BoundaryOutput;
						break;
			default:	fileName = "";
		}
		System.out.println("filename: "+fileName);
		FileWriter fw = new FileWriter(fileName);
		BufferedWriter bw = new BufferedWriter(fw);
		//write header
		System.out.println("Header: "+boundaryHeader);
		bw.write(boundaryHeader);
		bw.newLine();
		
		int count = 0;
		for(double radius : radiiForBoundary){
			System.out.println(count++);
			//reset neighborhood values for significance testing
			neighborCount=0;
			double marker1IF = getAveIfForRadius(MARKER1, ifType, radius);
			double marker2IF = getAveIfForRadius(MARKER2, ifType, radius);
			double percDiff = (marker1IF-marker2IF)/marker2IF*100;
			double pVal = calculatePValue(ifType);
			double avgNeighborhoodCells = ((double)neighborCount)/(popIndCellsMarker1.size()+popIndCellsMarker2.size());
			System.out.println("num of neighborhood cells: "+avgNeighborhoodCells);
			bw.write(radius+","+marker1IF+","+marker2IF+","+percDiff+","+pVal+","+avgNeighborhoodCells);
			bw.newLine();
		}
		bw.close();
		fw.close();
		System.out.println(fileName+" complete!");
	}
	
	
	private static void writePopulationComparisonFile() throws IOException{
		FileWriter fw = new FileWriter(populationOutput);
		BufferedWriter bw = new BufferedWriter(fw);
		//write header
		bw.write(populationHeader);
		bw.newLine();
		//cycle through each population sample
		for(int i=0; i<indMeans.size(); i++){
			getPopIndCells(indMeanLows.get(i), indMeanHighs.get(i));
			bw.write(indMeans.get(i)+","+popIndCellsMarker1.size()+","+popIndCellsMarker2.size());
			bw.newLine();
		}
		bw.close();
		fw.close();
		System.out.println(populationOutput+" complete!");
	}
	
	
	private static void writeBoundarySummaryFile() throws IOException{
		FileWriter fw = new FileWriter(boundarySummaryOutput);
		BufferedWriter bw = new BufferedWriter(fw);
		
		//read csvs to find max differences and their radii
		BufferedReader br1 = new BufferedReader(new FileReader(if1BoundaryOutput));
		double max1=-10;
		double radius1=0;
		String line;
		int lineNum = 0;
		while((line = br1.readLine()) != null){
			//skip the header
			if(lineNum>0){
				String[] vals = line.split(",");
				if(Double.parseDouble(vals[3])>max1){
					max1 = Double.parseDouble(vals[3]);
					radius1 = Double.parseDouble(vals[0]);
				}
			}
			lineNum++;
		}
		br1.close();
		
		BufferedReader br2 = new BufferedReader(new FileReader(if2BoundaryOutput));
		double max2=-10;
		double radius2=0;
		lineNum = 0;
		while((line = br2.readLine()) != null){
			//skip the header
			if(lineNum>0){
				String[] vals = line.split(",");
				if(Double.parseDouble(vals[3])>max2){
					max2 = Double.parseDouble(vals[3]);
					radius2 = Double.parseDouble(vals[0]);
				}
			}
			lineNum++;
		}
		br2.close();
		
		BufferedReader br3 = new BufferedReader(new FileReader(if3BoundaryOutput));
		double max3=-10;
		double radius3=0;
		lineNum=0;
		while((line = br3.readLine()) != null){
			//skip the header
			if(lineNum>0){
				String[] vals = line.split(",");
				if(Double.parseDouble(vals[3])>max3){
					max3 = Double.parseDouble(vals[3]);
					radius3 = Double.parseDouble(vals[0]);
				}
			}
			lineNum++;
		}
		br3.close();
		
		BufferedReader br4 = new BufferedReader(new FileReader(if4BoundaryOutput));
		double max4=-10;
		double radius4=0;
		lineNum=0;
		while((line = br4.readLine()) != null){
			//skip the header
			if(lineNum>0){
				String[] vals = line.split(",");
//				System.out.println(vals[3]);
				if(Double.parseDouble(vals[3])>max4){
					max4 = Double.parseDouble(vals[3]);
					radius4 = Double.parseDouble(vals[0]);
				}
			}
			lineNum++;
		}
		br4.close();
		
		//write header
		bw.write(boundarySummaryHeader);
		bw.newLine();
		//if 1
		bw.write(IF1+", "+max1+", "+radius1);
		bw.newLine();
		//if 2
		bw.write(IF2+", "+max2+", "+radius2);
		bw.newLine();
		//if 1
		bw.write(IF3+", "+max3+", "+radius3);
		bw.newLine();
		//if 1
		bw.write(IF4+", "+max4+", "+radius4);
		
		bw.close();
		fw.close();
		System.out.println(boundarySummaryOutput+" complete!");
	}
	
	
	private static Connection getConnection() throws SQLException{
		if(c == null){
			try {
				Class.forName("org.sqlite.JDBC");
			} catch (ClassNotFoundException e) {
				System.err.println("Could not find Driver: org.sqlite.JDBC");
				e.printStackTrace();
			}
			c = DriverManager.getConnection(dbURL);
		}
		return c;
	}
	
	private static void cleanUp(Statement s, ResultSet rs){
		try {
			if(s!=null){
				s.close();
			}
			if(rs!=null){
				rs.close();
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}
	
	private static ArrayList<Cell> getAllCells(){
		ArrayList<Cell> cells = new ArrayList<Cell>();
		
		Statement s = null;
		ResultSet rs = null;
		
		try{
			c = getConnection();
			s = c.createStatement();
			
//			query = "select count(*) from per_object";
			query = "select objectNumber, imageNumber, nuclei_location_center_x, "
					+ "nuclei_location_center_y, nuclei_math_normalized_"+indType+", "
					+ "nuclei_math_normalized_"+MARKER1+", nuclei_math_normalized_"+MARKER2
					+ " from per_object";
			
			rs = s.executeQuery(query);
			
			while(rs.next()){
				int id = rs.getInt("objectNumber");
				int imageId = rs.getInt("imageNumber");
				double x = rs.getDouble("nuclei_location_center_x");
				double y = rs.getDouble("nuclei_location_center_y");
				double ind = rs.getDouble("nuclei_math_normalized_"+indType);
				double m1 = rs.getDouble("nuclei_math_normalized_"+MARKER1);
				double m2 = rs.getDouble("nuclei_math_normalized_"+MARKER2);
				
				cells.add(new Cell(id, imageId, x, y, ind, m1, m2));			}
			
		}catch(SQLException e){
			System.err.println(query);
			System.err.println(e);
		}finally{
			cleanUp(s, rs);
		}
		
		return cells;
	}
	
	private static int getNumberOfCellType(String cellType, String cellIntensity){
		int result = -1;
		
		Statement s = null;
		ResultSet rs = null;
		
		try{
			c = getConnection();
			s = c.createStatement();
			
			query = "select count(*) from per_object where "+cellType+" > 0 and "+cellIntensity+">0";
			
			rs = s.executeQuery(query);
			
			while(rs.next()){
				result = rs.getInt("count(*)");
			}
			
		}catch(SQLException e){
			
		}finally{
			cleanUp(s, rs);
		}
		
		return result;
	}
	
	private static double getMeanInd(){
		double result = -1;
		
		Statement s = null;
		ResultSet rs = null;
		
		try{
			c = getConnection();
			s = c.createStatement();
			
			query = "select avg("+indIntensity+"), max("+indIntensity+") from per_object where "+indCellType+" > 0 and "+indIntensity+">0";
			
			rs = s.executeQuery(query);
			
			while(rs.next()){
				result = rs.getDouble("avg("+indIntensity+")");
				indMax = rs.getDouble("max("+indIntensity+")");
			}
		}catch(SQLException e){
			
		}finally{
			cleanUp(s, rs);
		}
		return result;
	}
	
	private static ArrayList<Cell> getPopIndCells(double lowBoundary, double highBoundary){
		ArrayList<Cell> indCells = new ArrayList<Cell>();
		
		//clear population arrays at beginning of each run
		popIndCellsBoth.clear();
		popIndCellsMarker1.clear();
		popIndCellsMarker2.clear();
		popIndCellsNone.clear();
		
//		System.out.println("cd34 size: "+meanHaCellsCd34.size());
		
		Statement s = null;
		ResultSet rs = null;
		
		try{
			c = getConnection();
			s = c.createStatement();
			
			String pre = "nuclei_math_normalized_";
			String indStr = pre+indType;
			String m1Str = pre+MARKER1;
			String m2Str = pre+MARKER2;
			
			
			query = "select objectNumber, imageNumber, nuclei_location_center_x, "
					+ "nuclei_location_center_y, "+indStr+", "
					+ m1Str+", "+m2Str+" "+"from per_object where "
					+ indIntensity+ ">"+lowBoundary+" and "+indIntensity+"<"+highBoundary;
			
//			System.out.println(query);
			
			rs = s.executeQuery(query);
			
			while(rs.next()){
				int id = rs.getInt("objectNumber");
				int i_id = rs.getInt("imageNumber");
				double x = rs.getDouble("nuclei_location_center_x");
				double y = rs.getDouble("nuclei_location_center_y");
				double ind = rs.getDouble(indStr);
				double m1 = rs.getDouble(m1Str);
				double m2 = rs.getDouble(m2Str);
				
				Cell cell = new Cell(id, i_id, x, y, ind, m1, m2);
				
				indCells.add(cell);
				
				//also sort into cebpa, cd34, and both
				if(cell.isMarker1Pos() && cell.isMarker2Pos()){
					popIndCellsBoth.add(cell);
				}else if(cell.isMarker1Pos()){
					popIndCellsMarker1.add(cell);
				}else if(cell.isMarker2Pos()){
					popIndCellsMarker2.add(cell);
				}else{
					popIndCellsNone.add(cell);
				}
			}
			
		}catch(SQLException e){
			
		}finally{
			cleanUp(s, rs);
		}
		
		return indCells;
	}
	
	
	private static void setIFOfTargetCell(String ifType, Cell targetCell, double radius){
		double ifValue = 0;
		
		int count = 0;
		//calculate distance from target cell to gather neighborhood cells
		double targetX = targetCell.getX();
		double targetY = targetCell.getY();
		int targetImage = targetCell.getImageId();
		for(Cell c : allCells){
			//if it is the same cell as the target cell, continue
			//if not from the same image, continue
			if(c.getID() == targetCell.getID() || c.getImageId() != targetImage){
				continue;
			}
			double x = c.getX();
			double y = c.getY();
			
			double xDiff = Math.abs(x-targetX);
			double yDiff = Math.abs(y-targetY);
			double dist = Math.sqrt(xDiff*xDiff+yDiff*yDiff);

			if(dist<radius){
				//do something different depending on ifType
				if(ifType.equals(IF1) || ifType.equals(IF2)){
					ifValue += c.getIndependentStainExpression();
					count++;
				}
				if(ifType.equals(IF3)){
					ifValue += dist*c.getIndependentStainExpression();
				}
				if(ifType.equals(IF4) && dist>0){
					ifValue += c.getIndependentStainExpression()/dist;
				}
				
				//increment the neighbor counter
				neighborCount++;
				//add to the marker neighborhood array for significance
			}
		}
		//we need to divide the total expression by cell count for if2
		if(ifType.equals(IF2)){
			ifValue = ifValue/count;
		}
		
		targetCell.setIFValue(ifValue, ifType);
	}

	
	/**
	 * Cycles through all CEBPa+ mean HA cells and calculates the IF
	 * for each one for a given radius.  Then returns the mean of that
	 * IF.
	 * @param radius
	 * @return
	 */
	private static double getAveIfForRadius(String cellGroup, String ifVersion, double radius){
		double result = 0;
		
//		ArrayList<Double> impactFactors = new ArrayList<Double>();
		
		ArrayList<Cell> cellsToInvestigate = new ArrayList<Cell>();
		if(cellGroup.equals(MARKER1)){
			cellsToInvestigate = new ArrayList<Cell>(popIndCellsMarker1);
		}
		else if(cellGroup.equals(MARKER2)){
			cellsToInvestigate = new ArrayList<Cell>(popIndCellsMarker2);
		}
		
		//get all the IFs for the mean ha cebpa cells
		for(Cell cell : cellsToInvestigate){
			setIFOfTargetCell(ifVersion, cell, radius);
		}
		
		//find the average
		double totalIF = 0;
		for(Cell c : cellsToInvestigate){
			totalIF += c.getIF(ifVersion);
		}
		
		result = totalIF/cellsToInvestigate.size();
		
		return result;
	}
}