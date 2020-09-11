import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;

import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RVector;
import org.rosuda.JRI.Rengine;

public class Main {
	private static String dbPath;
	private static String dbURL;
	private static Connection c;
	private static String query;
	
	private static String outputPath;
	private static String cell1ArrayOutputPath;
	private static String cell2ArrayOutputPath;
	
	private static String cell1Name = "ckHA"; //cd34
	private static String cell2Name = "msHA"; //cebpa
	//The cell type that is driving the change (HA, GATA6, etc.)
	private static String cell3Name = "GATA6";
	
	private static String cell3Intensity = "nuclei_math_normalized_"+cell3Name;
	private static String cell1Intensity = "nuclei_math_normalized_"+cell1Name;
	private static String cell2Intensity = "nuclei_math_normalized_"+cell2Name;
	private static String cell3Type = "nuclei_children_"+cell3Name+"_stain_count";
	private static String cell1Type = "nuclei_children_"+cell1Name+"_stain_count";
	private static String cell2Type = "nuclei_children_"+cell2Name+"_stain_count";
	private static String meanStr = "avg";
	private static String medianStr = "median";
	
	private static int numOfCells;
	private static int numOfCell3Cells;
	private static int numOfCell1Cells;
	private static int numOfCell2Cells;
	
	private static double meanCell3;
	private static double meanCell1;
	private static double meanCell2;
	
	private static double medianCell3;
	private static double medianCell1;
	private static double medianCell2;
	
	private static int numOfCell1CellsLessThanMeanCell3;
	private static int numOfCell2CellsLessThanMeanCell3;
	private static int numOfCell1CellsGreaterThanMeanCell3;
	private static int numOfCell2CellsGreaterThanMeanCell3;
	private static int numOfCell1CellsWithNoCell3;
	private static int numOfCell2CellsWithNoCell3;
	
	private static int numOfCell2CellsWithNoCell1Cells;
	
	private static double meanCell3ForCell1Cells;
	private static double meanCell3ForCell2Cells;
	private static double medianCell3ForCell1Cells;
	private static double medianCell3ForCell2Cells;
	
	private static String decFormat = "%.2f";
	
	public static void main(String[] args){
		
//		String basePath = "C:/Users/shay_/Documents/ASU/thesis/Images for Analysis/Early Time Points/Day 3/PGP1 D3, T, GATA6, FOXA2, reimaged 11-11-16/cropped/standard_crop2/output/";
//		String basePath = "C:/Users/shay_/Documents/ASU/thesis/Images for Analysis/thesis stainings/52 - PGP1 D6 - msSOX2, rbT, gtGATA6/cropped/standard_crop2/output/";
//		String basePath = "C:/Users/shay_/Documents/ASU/thesis/Images for Analysis/Early Time Points/Day 3/PGP1 D3, T, GATA6, FOXA2, reimaged 11-11-16/cropped/standard_crop2/output/";
//		String basePath = "C:/Users/shay_/Documents/ASU/thesis/Images for Analysis/CEBPa-CD34-HA-DAPI/cropped/standard_crop2/output/";
		String basePath = "C:/Users/shay_/Dropbox (ASU)/GATA6 Expression Analysis/Cell Profiler Database Files/55 - PGP1-HA D3, ckHAg, msHAr, rbGATA6fr/";
		
		dbPath = basePath+"55 - PGP1-HA D3, ckHAg, msHAr, rbGATA6fr, subtract 1_1.db";
		System.out.println(dbPath);
//		dbPath = basePath+"standard_crop2.db";
		outputPath = basePath+"pgp1-ha-d3-ckHA-msHA-gata6-subtract-1_1-results.txt";
//		outputPath = basePath+"pgp1-d14-cd34-cebpa-ha-results.txt";
		
		cell1ArrayOutputPath = basePath+cell1Name+"_"+cell3Name+"_values.txt";
		cell2ArrayOutputPath = basePath+cell2Name+"_"+cell3Name+"_values.txt";

		dbURL = "jdbc:sqlite:"+dbPath;
		
		
		//populate stats
		numOfCells = getTotalNumberOfCells();
		numOfCell3Cells = getNumberOfCellType(cell3Type);
		numOfCell1Cells = getNumberOfCellType(cell1Type);
		numOfCell2Cells = getNumberOfCellType(cell2Type);
		meanCell3 = getMeanOrMedianValue(meanStr, cell3Intensity);
		meanCell1 = getMeanOrMedianValue(meanStr, cell1Intensity);
		meanCell2 = getMeanOrMedianValue(meanStr, cell2Intensity);
		medianCell3 = getMeanOrMedianValue(medianStr, cell3Intensity);
		medianCell1 = getMeanOrMedianValue(medianStr, cell1Intensity);
		medianCell2 = getMeanOrMedianValue(medianStr, cell2Intensity);
		numOfCell1CellsGreaterThanMeanCell3 = getCellTypeAboveOrBelowCell3Mean(cell1Type, ">");
		numOfCell2CellsGreaterThanMeanCell3 = getCellTypeAboveOrBelowCell3Mean(cell2Type, ">");
		numOfCell1CellsLessThanMeanCell3 = getCellTypeAboveOrBelowCell3Mean(cell1Type, "<");
		numOfCell2CellsLessThanMeanCell3 = getCellTypeAboveOrBelowCell3Mean(cell2Type, "<");
		numOfCell1CellsWithNoCell3 = getCellNumberNoCell3(cell1Type);
		numOfCell2CellsWithNoCell3 = getCellNumberNoCell3(cell2Type);
		
		numOfCell2CellsWithNoCell1Cells = getCellNumberNoCellType(cell2Type, cell1Type);
		
		meanCell3ForCell1Cells = getMeanOrMedianCell3ValueForCellType(cell1Type, meanStr);
		meanCell3ForCell2Cells = getMeanOrMedianCell3ValueForCellType(cell2Type, meanStr);
		
		medianCell3ForCell1Cells = getMeanOrMedianCell3ValueForCellType(cell1Type, medianStr);
		medianCell3ForCell2Cells = getMeanOrMedianCell3ValueForCellType(cell2Type, medianStr);
		
		//TODO calculate geometric mean for cell populations: exp(average(ln(data))
		
		//calculate p value for mean cell3 type (gata6) difference between cell1 and cell2
//		double pValue = calculatePValue();
//		System.out.println(pValue);
		
		//write out to a file
		try {
			FileWriter fw = new FileWriter(outputPath);
			BufferedWriter bw = new BufferedWriter(fw);
			
			System.out.println("Writing summary report...");
			generateReport(bw);
		
			bw.close();
			fw.close();
//			//write out the array files to use in r for the t-tests
//			FileWriter fw1 = new FileWriter(cell1ArrayOutputPath);
//			BufferedWriter bw1 = new BufferedWriter(fw1);
//			System.out.println("Writing "+cell1Name+" array file...");
//			generateArrayTextFiles(bw1, cell1Type);
//			
//			FileWriter fw2 = new FileWriter(cell2ArrayOutputPath);
//			BufferedWriter bw2 = new BufferedWriter(fw2);
//			System.out.println("Writing "+cell2Name+" array file...");
//			generateArrayTextFiles(bw2, cell2Type);
//			
//			bw1.close();
//			fw1.close();
//			bw2.close();
//			fw2.close();
			
			
		} catch (IOException e1) {
			System.err.println("Error writing file.");
			e1.printStackTrace();
		}
		
		
		
		//close connection
		try {
			c.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
		
		System.out.println("All files generated successfully!");
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
	
	private static int getTotalNumberOfCells(){
		int result = -1;
		
		Statement s = null;
		ResultSet rs = null;
		
		try{
			c = getConnection();
			s = c.createStatement();
			
			query = "select count(*) from per_object";
			
			rs = s.executeQuery(query);
			
			while(rs.next()){
				result = rs.getInt("count(*)");
			}
			
		}catch(SQLException e){
			System.err.println(query);
			System.err.println(e);
		}finally{
			cleanUp(s, rs);
		}
		
		return result;
	}
	
	private static int getNumberOfCellType(String cellType){
		int result = -1;
		
		Statement s = null;
		ResultSet rs = null;
		
		try{
			c = getConnection();
			s = c.createStatement();
			
			query = "select count(*) from per_object where "+cellType+" > 0";
//			System.out.println(query);
			
			rs = s.executeQuery(query);
			
			
			while(rs.next()){
//				System.out.println("here");
				result = rs.getInt("count(*)");
//				System.out.println(result);
			}
			
		}catch(SQLException e){
			
		}finally{
			cleanUp(s, rs);
		}
		
		return result;
	}
	
	private static double getMeanOrMedianValue(String meanOrMedian, String intensityType){
		double result = -1;
		
		Statement s = null;
		ResultSet rs = null;
		
		try{
			c = getConnection();
			s = c.createStatement();
			
			query = "select "+meanOrMedian+"("+intensityType+") from per_object;";
//			System.out.println(query);
			
			rs = s.executeQuery(query);
			
			while(rs.next()){
				result = rs.getDouble(meanOrMedian+"("+intensityType+")");
//				System.out.println(result);
			}
			
		}catch(SQLException e){
			System.err.println(query);
			e.printStackTrace();
		}finally{
			cleanUp(s, rs);
		}
		
		return result;
	}
	
	
	private static int getCellTypeAboveOrBelowCell3Mean(String cellType, String aboveBelow){
		int result = -1;
		
		Statement s = null;
		ResultSet rs = null;
		
		try{
			c = getConnection();
			s = c.createStatement();
			
			query = "select count(*) from per_object where "+cellType+"> 0 and "
					+ cell3Intensity +" "+aboveBelow+" "+meanCell3;
			
//			System.out.println(query);
			rs = s.executeQuery(query);
			
			while(rs.next()){
				result = rs.getInt("count(*)");
			}
		}catch(SQLException e){
			System.err.println(query);
			e.printStackTrace();
		}finally{
			cleanUp(s, rs);
		}
		
		return result;
	}
	
	
	private static int getCellNumberNoCell3(String cellType){
		int result = -1;
		
		Statement s = null;
		ResultSet rs = null;
		
		try{
			c = getConnection();
			s = c.createStatement();
			
			query = "select count(*) from per_object where "+cellType+"> 0 and "
					+ cell3Type + "=0";
			
//			System.out.println(query);
			rs = s.executeQuery(query);
			
			while(rs.next()){
				result = rs.getInt("count(*)");
//				System.out.println(result);
			}
		}catch(SQLException e){
			System.err.println(query);
			e.printStackTrace();
		}finally{
			cleanUp(s, rs);
		}
		
		return result;
	}
	
	
	private static int getCellNumberNoCellType(String numberCellType, String noCellType){
		int result = -1;
		
		Statement s = null;
		ResultSet rs = null;
		
		try{
			c = getConnection();
			s = c.createStatement();
			
			query = "select count(*) from per_object where "+numberCellType+"> 0 and "
					+ noCellType + "=0";
			
//			System.out.println(query);
			rs = s.executeQuery(query);
			
			while(rs.next()){
				result = rs.getInt("count(*)");
//				System.out.println(result);
			}
		}catch(SQLException e){
			System.err.println(query);
			e.printStackTrace();
		}finally{
			cleanUp(s, rs);
		}
		
		return result;
	}
	
	
	private static double getMeanOrMedianCell3ValueForCellType(String cellType, String meanOrMed){
		double result = -1;
		
		Statement s = null;
		ResultSet rs = null;
		
		try{
			c = getConnection();
			s = c.createStatement();
			
			//if avg...make sure to include cells that have 0 cell3Intensity
			query = "select ("+meanOrMed+"("+cell3Intensity+")*count("+cell3Intensity+")/count(objectNumber)) "
					+ "from per_object where "+cellType+">0";
//			System.out.println(query);
			
			rs = s.executeQuery(query);
			
			while(rs.next()){
				result = rs.getDouble(1);
//				System.out.println(result);
			}
			
		}catch(SQLException e){
			System.err.println(query);
			e.printStackTrace();
		}finally{
			cleanUp(s, rs);
		}
		
		return result;
	}
	
	
	private static double[] getAllCell3LevelsForCellType(String cellType){
		double[] result = null;
		int size = 0;
		
		Statement s = null;
		ResultSet rs = null;
		
		try{
			c = getConnection();
			s = c.createStatement();
			
			String query1 = "select count(objectNumber) from per_object where "+cellType+">0";
			ResultSet rs1 = s.executeQuery(query1);
			while(rs1.next()){
				size = rs1.getInt(1);
			}
			
			result = new double[size];
			
			query = "select objectNumber, "+cell3Intensity+" from per_object "
					+ "where "+cellType+">0";
			
			rs = s.executeQuery(query);
			
			int i = 0;
			while(rs.next()){
				result[i] = rs.getDouble(2);
				i++;
			}
			
		}catch(SQLException e){
			System.err.println(query);
			e.printStackTrace();
		}finally{
			cleanUp(s, rs);
		}
		
		
		return result;
	}
	
	
	private static double calculatePValue(){
		double result =0;
		
		System.out.println("setting up r engine");
		Rengine engine = new Rengine(new String[] {"--no-save"}, false, null);
		
		String a1 = cell1Name+"_"+cell3Name;
		String a2 = cell2Name+"_"+cell3Name;
		
		engine.assign(a1, getAllCell3LevelsForCellType(cell1Type));
		engine.assign(a2, getAllCell3LevelsForCellType(cell2Type));
		
		REXP rexp = engine.eval("t.test("+a2+","+a1+", paired=FALSE)");
//		System.out.println(rexp.getContent());
		
		RVector rvec = rexp.asVector();
		String pString = rvec.get(2).toString();
		
		//Take off beginning: "[REAL* (" and end: ")]"
		pString = pString.substring(8, pString.length()-2);
		
		
//		System.out.println(pString);
		
		result = Double.parseDouble(pString);
		
//		System.out.println(engine.eval("t.test("+a2+","+a1+", paired=FALSE)"));
		
		engine.end();
		
		return result;
	}
	
	
	private static void generateReport(BufferedWriter bw) throws IOException{
		//title
		//TODO: let this be a variable the user can specify ahead of time
		bw.write(cell1Name+"-"+cell2Name+"-"+cell3Name+" Normalized Image Results Summary");
		bw.newLine();
		bw.newLine();
		//totals section
		bw.write("Number of identified cells (nuclei): "+numOfCells);
		bw.newLine();
		bw.write("Number of identified "+cell3Name+" cells: "+numOfCell3Cells);
		bw.newLine();
		bw.write("Number of identified "+cell1Name+" cells: "+numOfCell1Cells);
		bw.newLine();
		bw.write("Number of identified "+cell2Name+" cells: "+numOfCell2Cells);
		bw.newLine();
		bw.newLine();
		//percent of cells
		double percCellsCell3 = (double)numOfCell3Cells/(double)numOfCells*100;
		double percCellsCell1 = (double)numOfCell1Cells/(double)numOfCells*100;
		double percCellsCell2 = (double)numOfCell2Cells/(double)numOfCells*100;
		bw.write("Percent of cells that are "+cell3Name+"+: "+String.format(decFormat, percCellsCell3)+"%");
		bw.newLine();
		bw.write("Percent of cells that are "+cell1Name+"+: "+String.format(decFormat, percCellsCell1)+"%");
		bw.newLine();
		bw.write("Percent of cells that are "+cell2Name+"+: "+String.format(decFormat, percCellsCell2)+"%");
		bw.newLine();
		bw.newLine();
		//mean values
		bw.write("Mean "+cell3Name+" value: "+meanCell3);
		bw.newLine();
		bw.write("Mean "+cell1Name+" value: "+meanCell1);
		bw.newLine();
		bw.write("Mean "+cell2Name+" value: "+meanCell2);
		bw.newLine();
		bw.newLine();
		//median values
		bw.write("Median "+cell3Name+" value: "+medianCell3);
		bw.newLine();
		bw.write("Median +"+cell1Name+" value: "+medianCell1);
		bw.newLine();
		bw.write("Median "+cell2Name+" value: "+medianCell2);
		bw.newLine();
		bw.newLine();
		//child cells with ha greater than the mean
		double cell1PercOverCell3Mean = (double)numOfCell1CellsGreaterThanMeanCell3/(double)numOfCell1Cells*100;
		double cell2PercOverCell3Mean = (double)numOfCell2CellsGreaterThanMeanCell3/(double)numOfCell2Cells*100;
		bw.write(cell1Name+"+ cells with "+cell3Name+" > mean: "+numOfCell1CellsGreaterThanMeanCell3+" (number), "+String.format(decFormat, cell1PercOverCell3Mean)+"%");
		bw.newLine();
		bw.write(cell2Name+"+ cells with "+cell3Name+" > mean: "+numOfCell2CellsGreaterThanMeanCell3+" (number), "+String.format(decFormat, cell2PercOverCell3Mean)+"%");
		bw.newLine();
		bw.newLine();
		//child cells with ha less than the mean
		double cell1PercUnderCell3Mean = (double)numOfCell1CellsLessThanMeanCell3/(double)numOfCell1Cells*100;
		double cell2PercUnderCell3Mean = (double)numOfCell2CellsLessThanMeanCell3/(double)numOfCell2Cells*100;
		bw.write(cell1Name+"+ cells with "+cell3Name+" < mean: "+numOfCell1CellsLessThanMeanCell3+" (number), "+String.format(decFormat, cell1PercUnderCell3Mean)+"%");
		bw.newLine();
		bw.write(cell2Name+"+ cells with "+cell3Name+" < mean: "+numOfCell2CellsLessThanMeanCell3+" (number), "+String.format(decFormat, cell2PercUnderCell3Mean)+"%");
		bw.newLine();
		bw.newLine();
		//child cells with no ha
		double cell1PercNoCell3 = (double)numOfCell1CellsWithNoCell3/(double)numOfCell1Cells*100;
		double cell2PercNoCell3 = (double)numOfCell2CellsWithNoCell3/(double)numOfCell2Cells*100;
		bw.write(cell1Name+"+ cells with no "+cell3Name+": "+numOfCell1CellsWithNoCell3+" (number), "+String.format(decFormat, cell1PercNoCell3)+"%");
		bw.newLine();
		bw.write(cell2Name+"+ cells with no "+cell3Name+": "+numOfCell2CellsWithNoCell3+" (number), "+String.format(decFormat, cell2PercNoCell3)+"%");
		bw.newLine();
		bw.newLine();
		//child cell with no other child
		bw.write(cell2Name+"+ cells with no "+cell1Name+": "+numOfCell2CellsWithNoCell1Cells);
		bw.newLine();
		bw.newLine();
		//mean and median ha value for child cells
		bw.write("Mean "+cell3Name+" value for "+cell1Name+"+ cells: "+meanCell3ForCell1Cells);
		bw.newLine();
		bw.write("Mean "+cell3Name+" value for "+cell2Name+"+ cells: "+meanCell3ForCell2Cells);
		bw.newLine();
		bw.write("P value for mean differences: "+calculatePValue());
		bw.newLine();
		bw.newLine();
		bw.write("Median "+cell3Name+" value for "+cell1Name+"+ cells: "+medianCell3ForCell1Cells);
		bw.newLine();
		bw.write("Median "+cell3Name+" value for "+cell2Name+"+ cells: "+medianCell3ForCell2Cells);
		bw.newLine();
	}
	
//	private static void generateArrayTextFiles(BufferedWriter bw, String cellType) throws IOException{
//		for(Double val : getAllCell3LevelsForCellType(cellType)){
//			bw.write(val+"");
//			bw.newLine();
//		}
//	}
}
