import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Vector;

public class DataAccess {
	private static String dbPath;
	private static String dbURL;
	private static Connection c;
	private static String query;
	
	private static String tableName = "per_object";
	
	
	public static void setDatabasePath(String path){
		dbPath = path;
		dbURL = "jdbc:sqlite:"+dbPath;
	}
	
	private static Connection getConnection() throws SQLException{
		if(c == null && dbURL != null){
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
	
	
	public static Vector<String> getColumnNames(){
		Vector<String> result = new Vector<String>();
		
		Statement s = null;
		ResultSet rs = null;
		
		try{
			c = getConnection();
			s = c.createStatement();
			
			query = "PRAGMA table_info("+tableName+")";
			
			rs = s.executeQuery(query);
			
			while(rs.next()){
				result.add(rs.getString("name"));
			}
		}catch(SQLException e){
			System.err.println(query);
			e.printStackTrace();
		}catch(NullPointerException e1){
			
		}finally{
			cleanUp(s, rs);
		}
		
		return result;
	}
}
