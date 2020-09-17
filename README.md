# Neighborhood-Impact-Factor

The concept of Impact Factor in the context of cell signaling attempts to integrate overlapping layers of local signals that influence a cell’s GRN into a simplified quantitative form. The radius at which a local signaling element influences the cell fate in the microenvironment is referred to here as a cell’s “neighborhood”. Neighborhood Impact Factor (NIF) is the quantitative assessment of a given signaling element within a certain neighborhood radius around a target cell.

We have proposed three ways to define the NIF, which we have termed “Total Expression”, “Local Density”, and “Distance Adjusted”. However, the models of neighborhood expression should be tailored for the analysis of the system of interest, depending on a proposed or expected expression pattern.

- **Total Expression Impact Factor**

  The impact factor of the Total Expression method is defined as:

  This method simply calculates the sum of the expression levels of all cells expressing a Gene of Interest (GOI) within the neighborhood by finding all the cells within the neighborhood radius and adding up their GOI expression levels as quantified by intensity. This expression level is then normalized to the gene expression of the cell for which the neighborhood is being analyzed.

- **Local Density Impact Factor**

  The Local Density method defines impact factor as:

  This method takes the sum of the expression of a given GOI within the neighborhood divided by all the cells in the neighborhood to produce an average expression per cell. This method can be useful to check if a density of nearby cells expressing a GOI, not just a threshold of local expression, is a key determining factor in driving the biological events associated with the cell neighborhood.

- **Distance Adjusted Impact Factor**

  The Distance Adjusted method defines impact factor as:

  This method quantifies the neighborhood with the assumption that the neighborhood cells expressing a GOI farther away from the cell of interest have reduced influence on that cell. Therefore, the distance away from the central cell of interest is inversely proportional to the impact of that cell on the NIF.

The best NIF for representation of signaling effects within a given system should be iteratively optimized through modeling and experimentation to fit the experimental parameters of that system.

This GOI is also referred to as the "independent" gene and the other genes (which are used to determine what fate the cell has taken) are also referred to as the "dependent" genes.


## CellProfilerSummaryApp

This java program is written to work with a database (\*.db) file that has been outtputted from a CellProfiler pipeline.  
It uses the information in that database file to output three text files:

- **Summary File**
  
  This is a summary file that has all kinds of statistics about the populations being analyzed.  It contains popuation totals, percents, mean and median values, comparisons across the "dependent" populations, mean and median "independent" gene value for the two populations, and the p-value for those mean differences.
  
- **Raw Value Files**
 
  Two more files will be output -- one for each of the "dependent" genes.  These files will contain the raw "independent" gene levels for every "dependent" gene of that type.  
  > For instance, if HA is the "independent" gene and CEBP&alpha; and CD34 are the "dependent" genes, then there will be two output files: cebpa_HA_values.txt and cd34_HA_values.txt.  The first file will contain the HA levels for all positively identified CEBP&alpha; cells, and the second will contain the HA levels for all positively identified CD34 cells.

## GeneralNeighborhoodAnalysis
