import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Vector;

import javax.swing.AbstractAction;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.UIManager;
import javax.swing.filechooser.FileFilter;

public class MainFrame extends JFrame{
	private JButton inputBtn;
	private JFileChooser inputChooser;
	private JTextField inputTF;
	private JLabel nucleiLbl;
	private JComboBox<String> nucleiBx;
	private JLabel neighborhoodLbl;
	private JComboBox<String> neighborhoodBx;
	private JLabel fate1Lbl;
	private JComboBox<String> fate1Bx;
	private JLabel fate2Lbl;
	private JComboBox<String> fate2Bx;
	private JPanel boxPnl;
	private JButton statOutputBtn;
	private JFileChooser statOutputChooser;
	private JTextField statOutputTF;
	private JButton statBtn;
	
	private JTabbedPane mainPane;
	private JPanel setupPnl;
	private JPanel statPnl;
	private JPanel nifPnl;
	
	private Color lightBlue = UIManager.getColor("TabbedPane.selected");
	private int pad = 5;
	private Insets in = new Insets(pad, pad, pad, pad);
	
	public MainFrame(){
		//set title
		this.setTitle("Neighborhood Analysis");
		
		//create tabbed pane that houses everything
		mainPane = new JTabbedPane(JTabbedPane.BOTTOM);
		
		//create and add tabs
		setupPnl = buildSetupPanel();
		statPnl = buildStatPanel();
		mainPane.addTab("Setup", setupPnl);
		mainPane.addTab("Summary Stats", statPnl);
		
		//add mainPane to this frame
		this.setContentPane(mainPane);
		
		//set size
		Dimension size = new Dimension(650, 400);
		setPreferredSize(size);
		setMinimumSize(size);
		
		//set visible
		this.setVisible(true);
	}
	
	
	private JPanel buildSetupPanel(){
		//create components
		JPanel panel = new JPanel();
		panel.setBackground(lightBlue);
		panel.setLayout(new GridBagLayout());
		
		inputBtn = new JButton(inputAct);
		//limit the file type
		inputChooser = new JFileChooser();
		inputChooser.setCurrentDirectory(new File(System.getProperty("user.home")));
		inputChooser.setDialogTitle("Choose Input File");
		FileFilter dbFilter = new FileFilter() {
			public boolean accept(File pathname) {
				if(pathname.isDirectory()){
					return true;
				}else{
					String path = pathname.getAbsolutePath();
					if(path.endsWith(".db")){
						return true;
					}
				}
				return false;
			}

			public String getDescription() {
				return "Database File";
			}
		};
		inputChooser.setFileFilter(dbFilter);
		inputTF = new JTextField(50);
		inputTF.setEditable(false);
		
		nucleiLbl = new JLabel("Nuclei Marker:");
		neighborhoodLbl = new JLabel("Neighborhood Marker:");
		fate1Lbl = new JLabel("Cell Fate 1 Marker:");
		fate2Lbl = new JLabel("Cell Fate 2 Marker:");
		boxPnl = new JPanel(new GridBagLayout());
		boxPnl.setBackground(lightBlue);
		refreshStainBoxes();
		
		//build layout
		int row = 0;
		panel.add(inputBtn, new GridBagConstraints(0, row, 2, 1, 0, 0, GridBagConstraints.CENTER, GridBagConstraints.NONE, in, pad, pad));
		panel.add(inputTF, new GridBagConstraints(0, ++row, 2, 1, 0, 0, GridBagConstraints.CENTER, GridBagConstraints.NONE, in, pad, pad));
		panel.add(boxPnl, new GridBagConstraints(0, ++row, 2, 4, 0, 0, GridBagConstraints.CENTER, GridBagConstraints.NONE, in, pad, pad));
		
		return panel;
	}
	
	private AbstractAction inputAct = new AbstractAction("Set Source File:") {
		public void actionPerformed(ActionEvent e) {
			int val = inputChooser.showDialog(MainFrame.this, "Okay");
			if(val == JFileChooser.APPROVE_OPTION){
				//set the path
				String path = inputChooser.getSelectedFile().getPath();
				DataAccess.setDatabasePath(path);
				//update the text field
				inputTF.setText(path);
				//update the comboboxes.
				refreshStainBoxes();
			}
		}
	};
	
	private void refreshStainBoxes(){
		//remove everything from the box panel
		boxPnl.removeAll();
		//get the column names and populate comboboxes
		Vector<String> colNames = DataAccess.getColumnNames();
		if(colNames.size()>0){
			nucleiBx = new JComboBox<String>(colNames);
			nucleiBx.setEnabled(true);
			neighborhoodBx = new JComboBox<String>(colNames);
			neighborhoodBx.setEnabled(true);
			fate1Bx = new JComboBox<String>(colNames);
			fate1Bx.setEnabled(true);
			fate2Bx = new JComboBox<String>(colNames);
			fate2Bx.setEnabled(true);
		}else{
			Dimension boxSize = new Dimension(200, 25);
			nucleiBx = new JComboBox<String>();
			nucleiBx.setEnabled(false);
			nucleiBx.setPreferredSize(boxSize);
			neighborhoodBx = new JComboBox<String>();
			neighborhoodBx.setEnabled(false);
			neighborhoodBx.setPreferredSize(boxSize);
			fate1Bx = new JComboBox<String>();
			fate1Bx.setEnabled(false);
			fate1Bx.setPreferredSize(boxSize);
			fate2Bx = new JComboBox<String>();
			fate2Bx.setEnabled(false);
			fate2Bx.setPreferredSize(boxSize);
		}
		//add everything back to the panel
		int row = 0;
		boxPnl.add(nucleiLbl, new GridBagConstraints(0, row, 1, 1, 0, 0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, in, pad, pad));
		boxPnl.add(nucleiBx, new GridBagConstraints(1, row, 1, 1, 0, 0, GridBagConstraints.LINE_START, GridBagConstraints.NONE, in, pad, pad));
		boxPnl.add(neighborhoodLbl, new GridBagConstraints(0, ++row, 1, 1, 0, 0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, in, pad, pad));
		boxPnl.add(neighborhoodBx, new GridBagConstraints(1, row, 1, 1, 0, 0, GridBagConstraints.LINE_START, GridBagConstraints.NONE, in, pad, pad));
		boxPnl.add(fate1Lbl, new GridBagConstraints(0, ++row, 1, 1, 0, 0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, in, pad, pad));
		boxPnl.add(fate1Bx, new GridBagConstraints(1, row, 1, 1, 0, 0, GridBagConstraints.LINE_START, GridBagConstraints.NONE, in, pad, pad));
		boxPnl.add(fate2Lbl, new GridBagConstraints(0, ++row, 1, 1, 0, 0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, in, pad, pad));
		boxPnl.add(fate2Bx, new GridBagConstraints(1, row, 1, 1, 0, 0, GridBagConstraints.LINE_START, GridBagConstraints.NONE, in, pad, pad));
		
		//revalidate
		boxPnl.revalidate();
	}
	
	
	private JPanel buildStatPanel(){
		JPanel panel = new JPanel(new GridBagLayout());
		panel.setBackground(lightBlue);
	
		
		
		return panel;
	}
}
