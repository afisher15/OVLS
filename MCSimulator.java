/** 
 * A Monte Carlo simulator for crystal growth.
 * 
 * Ali Fisher 2015
 * 
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;

public class MCSimulator {
	/** A pixel map of the crystals. Every cell is a crystal id or 0 if empty space **/
	private short[][] crystalMap;         
	
	/** Counts for each pixel of cycles during which a monomer has been on the pixel **/
	private short[][] contourDiagram;
	
	/** A tessellation of the crystal space. Each entry is the id of the crystal to 
	 *  which the position is closest.
	 */
	private short[][] tessellation;
	
	/** A map from crystal id to crystal **/
	private Map<Short,Crystal> crystals;
	
	/** A list of monomers that have been used in this simulation **/
	private List<Monomer> monomers;
	
	/** Types of moves that a monomer can undergo, each with a designated probability **/
	public enum Move {
		WALK,          // move randomly in one of the 4 non-diagonal directions
		STICK_EDGE,    // stick to the edge of a crystal
		STICK_TOP,     // stick to the top of a crystal
		EDGE_BOUND,    // become edge-bound to a crystal
		UN_EDGE_BOUND; // become un-edge-bound from a crystal
		
		/**
		 * @return The probability of this move.
		 */
		double probability() {
		    switch(this) {
	            case WALK:   		return 1.0;
	            case STICK_EDGE:  	return 0.001;
	            case STICK_TOP:  	return 0;
	            case EDGE_BOUND: 	return 0;
	            case UN_EDGE_BOUND: return 0;
		    }
		    throw new AssertionError("Unknown move: " + this);
		}
	}
	
	// used to generate random numbers throughout the program
	private static final Random random = new Random();
	
	/**
	 * Creates a new a new MC simulator with the given crystal image.
	 * 
	 * @param filename A csv file containing 
	 */
	public MCSimulator(String filename) {
		crystalMap = arrayFromCSV(filename);
		//tessellation = computeTesselation(crystalMap);
		tessellation = arrayFromCSV("tessellation.csv");  // read in from file for efficiency
		contourDiagram = new short[crystalMap.length][crystalMap[0].length];
		crystals = initCrystals();
		monomers = new ArrayList<Monomer>();
	}
	
	/**
	 * Creates one monomer, places it randomly, and has it walk randomly until
	 * sticking to a crystal or falling off the edge of the image.
	 */
	public void doRandomWalk() {
		// create a new monomer and add it to the list of all monomers used in the
		// simulation
		Monomer m = new Monomer();
		monomers.add(m);
		
		// get the initial row and column for the monomer randomly
		int row = random.nextInt(crystalMap.length);
		int column = random.nextInt(crystalMap[0].length);
		
		// Figure out the crystal whose capture zone the monomer arrived in, and
		// increment this count for the crystal
		short captureZoneArrived = tessellation[row][column];
		crystals.get(captureZoneArrived).numArrived++;

		boolean stuck = false;     // whether on not the monomer has become stuck to a crystal
		boolean edgeBound = false; // whether or not the monomer is edge bound to a crystal
		boolean onMap = true;      // whether or not we have fallen off the map
		
		// The set of crystals that this monomer has been next to or on top of.
		// Used to keep track of sticking opportunities for crystals.
		Set<Short> crystalsTouched = new HashSet<Short>();
		
		// Walk until the monomer becomes stuck or walks off the edge
		while (!stuck && onMap) {
			// increment the contour diagram for this position
			contourDiagram[row][column]++;
			
			// get the possible moves for this monomer based on its position
			List<Move> possibleMoves = possibleMoves(row, column, edgeBound);
			
			// Add to our set of crystals that have had sticking opps during this MC cycle
			if (possibleMoves.contains(Move.STICK_TOP)) {
				crystalsTouched.add(crystalMap[row][column]);
			} else if (possibleMoves.contains(Move.STICK_EDGE)) {
				crystalsTouched.add(adjacentCrystal(row, column));
			}
			
			Move nextMove = chooseMove(possibleMoves, m);
			switch (nextMove) {
			 	case WALK:
			 		if (!edgeBound) {
				 		// choose randomly among the 4 directions to move
			 			int direction = random.nextInt(4);
				 		if (direction == 0) {
				 			row++;
				 		} else if (direction == 1) {
				 			column++;
				 		} else if (direction == 2) {
				 			row--;
				 		} else if (direction == 3) {
				 			column--;
				 		}
			 		} else {
			 			// TODO: fill in the edge bound move 
			 			// Calculate the two possible moves and decide randomly
			 			// between them
			 		}
			 		// increment the distance traveled for the monomer
			 		m.distT++;
			 		break;
	            case STICK_EDGE:
	            	stuck = true;
	            	break;
	            case STICK_TOP:
	            	stuck = true;
	            	break;
	            case EDGE_BOUND:
	            	edgeBound = true;
	            	break;
	            case UN_EDGE_BOUND:
	            	edgeBound = false;
	            	break;
	            default: 
	            	throw new RuntimeException("unknown move");	
			}
			
			/*if (row < 0 || row >= crystalMap.length || column < 0 || column >= crystalMap[0].length)
				// We walked off the edge! End this random walk
				onMap = false;*/
			if (row < 0)
				row = crystalMap.length - 1;
			else if (row > crystalMap.length - 1)
				row = 0;
			if (column < 0)
				column = crystalMap[0].length - 1;
			else if (column > crystalMap[0].length - 1)
				column = 0;
		}
		
		if (onMap) {
			// The monomer has stuck! Update this position in the crystal map
			// to reflect that it is now part of a crystal
			short adjacentCrystal = adjacentCrystal(row, column);
			
			// sanity check (since we're now stuck, we should be next to a crystal)
			assert adjacentCrystal != 1 : "no adjacent crystal found";
			
			// This spot now belongs to the adjacent crystal
			crystalMap[row][column] = adjacentCrystal;
			
			if (adjacentCrystal == captureZoneArrived)
				// we stuck to the same crystal whose capture zone we arrived in,
				// update this count for the crystal
				crystals.get(adjacentCrystal).numStuck++;
		}
		
		for (short cId: crystalsTouched) {
			// increment the # of sticking opps for any crystal this monomer touched
			crystals.get(cId).numStickingOps++;
		}
	}
	
	/**
	 * Returns a list of possible moves for the monomer based on its position and 
	 * whether or not it is edge bound.
	 * 
	 * @param row The monomer's row 
	 * @param column The monomer's column
	 * @param edgeBound True if the monomer is edge-bound
	 * @return A list of possible moves
	 */
	private List<Move> possibleMoves(int row, int column, boolean edgeBound) {
		List<Move> moves = new ArrayList<Move>();
		
		// It is always possible to walk so add this to the list
		moves.add(Move.WALK);
		
		if (crystalMap[row][column] != 0) {
			// we are on a crystal so it is possible to stick to the top
			// moves.add(Move.STICK_TOP);
		} else if (adjacentToCrystal(row, column)) {
			// we are next to a crystal; possible to stick to an edge
			moves.add(Move.STICK_EDGE);
			// also possible to become edge-bound if we aren't already
			/*if (!edgeBound)
				moves.add(Move.EDGE_BOUND);*/
		}
		
		// can become un-edge-bound if we are currently bound
		/*if (edgeBound)
			moves.add(Move.UN_EDGE_BOUND);*/
		
		return moves;
	}
	
	/**
	 * Pick a move for the random walker based on weighted probabilities of
	 * possible moves and update the monomer's time in existence.
	 * 
	 * @param possibleMoves The list of possible moves to choose from
	 * @param mmer The monomer to update the time for
	 * @return The chosen move
	 */
	private Move chooseMove(List<Move> possibleMoves, Monomer mmer) {
		// cumulative functions for possible moves
		List<Double> R_i = new ArrayList<Double>();
		double cumProb = 0.0;
		
		for (int i = 0; i < possibleMoves.size(); i++) {
			cumProb += possibleMoves.get(i).probability();
			R_i.add(cumProb);
		}
		
		double unif0to1 = random.nextDouble();
		double R_n = R_i.get(R_i.size() - 1);
		double uR_n = unif0to1 * R_n;
		
		// update the monomer's time based on the cumulative probability
		mmer.updateTime(R_n);
	
		for (int i = 0; i < possibleMoves.size(); i++) {
			double r_i = R_i.get(i);
			if (r_i >= uR_n)
				return possibleMoves.get(i);
		}
		
		return null;
	}
	
	/**
	 * Returns the id of the adjacent crystal if there is one. Returns -1 if the pixel
	 * at the given row/column does not have a face touching a crystal.
	 */
	private short adjacentCrystal(int row, int column) {
		// If we want to account for the possibility of multiple adjacent crystals,
		// this method may need to be changed later
		if (row != 0 && crystalMap[row - 1][column] != 0)
			return crystalMap[row - 1][column];
		if (row != crystalMap.length - 1 && crystalMap[row + 1][column] != 0)
			return crystalMap[row + 1][column];
		if (column != 0 && crystalMap[row][column - 1] != 0)
			return crystalMap[row][column - 1];
		if (column != crystalMap[0].length - 1 && crystalMap[row][column + 1] != 0)
			return crystalMap[row][column + 1];
	
		// no adjacent crystal found
		return -1;
	}
	
	/**
	 * Returns true if the pixel at the given row/column is adjacent to a crystal 
	 * (has at least one face touching a crystal pixel), false otherwise.
	 */
	private boolean adjacentToCrystal(int row, int column) {
		return adjacentCrystal(row, column) != -1;
	}
	
	/**
	 * Initializes the set of crystals in this image with their
	 * initial area, perimeter, and whether or not an edge crystal (edge crystal 
	 * being defined as one whose capture zone touches the end of the image)
	 * 
	 * @return A map of crystal ids to crystals
	 */
	private Map<Short,Crystal> initCrystals() {
		Map<Short,Crystal> crystals = new TreeMap<Short, Crystal>();
		for (int i = 0; i < crystalMap.length; i++) {
			for (int j = 0; j < crystalMap[0].length; j++) {
				short cId = crystalMap[i][j];
				if (cId == 0) 
					// no crystal here, advance to next pixel
					continue; 
				Crystal c;
				if (crystals.containsKey(cId)) {
					// we've already seen this crystal
					c = crystals.get(cId);
				} else {
					// create a new crystal
					c = new Crystal(cId);
					crystals.put(cId, c);
				}
				c.initArea++;  // increment the crystal area for this pixel
				//c.initPerim += numPerimEdges(i, j, cId);
				if (isEdgePixel(i, j, cId))
					// this pixel has a non-crystal pixel next to it; increment
					// the perimeter count
					c.initPerim++;				
			}
		}
	
		// Circle around the outer edge of the tessellation. Any crystal
		// IDs found along these edges correspond to edge crystals,
		// so mark them as such.
		for (int i = 0; i < tessellation.length; i++) {
			if (!crystals.containsKey(tessellation[i][0])) {
				System.out.println(tessellation[i][0]);
				System.out.println(i);
			}
			if (!crystals.containsKey(tessellation[i][tessellation[0].length - 1])) {
				System.out.println(tessellation[i][tessellation[0].length - 1]);
				System.out.println(i);
			}
			crystals.get(tessellation[i][0]).edge = true;
			crystals.get(tessellation[i][tessellation[0].length - 1]).edge = true;	
		}
		
		for (int i = 0; i < tessellation[0].length; i++) {
			crystals.get(tessellation[0][i]).edge = true;
			crystals.get(tessellation[tessellation.length - 1][i]).edge = true;	
		}		
		return crystals;
	}
	
	/**
	 * Computes the final areas and perimeters for the crystals. To be called
	 * at the end of the simulation before printing crystal stats.
	 */
	private void finalizeCrystals() {
		for (int i = 0; i < crystalMap.length; i++) {
			for (int j = 0; j < crystalMap[0].length; j++) {
				short cId = crystalMap[i][j];
				if (cId == 0) 
					// no crystal here, advance to next pixel
					continue; 
				Crystal c = crystals.get(cId);
				c.finalArea++;  // increment the crystal area for this pixel
				if (isEdgePixel(i, j, cId))
					// this pixel has a non-crystal pixel next to it; increment
					// the perimeter count
					c.finalPerim++;				
			}
		}

	}

	/**
	 * Returns true if this pixel has a face touching a pixel that is empty
	 * space, has an id is different from the given id,  empty space or the edge of the
	 * diagram.
	 * 
	 * @param row The row of the pixel
	 * @param column The column of the pixel
	 * @param cId The id of the crystal w
	 * @return true if an edge pixel
	 */
	private boolean isEdgePixel(int row, int column, short cId) {
		if (row == 0 || column == 0 || column == crystalMap[0].length - 1 ||
				row == crystalMap.length - 1) {
			// we are on the edge of the map!
			return true;
		}
		
		return  crystalMap[row][column - 1] != cId ||
				crystalMap[row][column + 1] != cId ||
				crystalMap[row - 1][column] != cId ||
				crystalMap[row - 1][column + 1] != cId ||
				crystalMap[row + 1][column] != cId ||
				crystalMap[row + 1][column - 1] != cId;
	}
	
	/////////////////////////////////////////////////////////////////////////////////////
	////                         RESULTS STATS METHODS                               ////
	/////////////////////////////////////////////////////////////////////////////////////
	
	/**
	 * Prints a CSV file of the image of the current simulation state.
	 * 
	 * @param filename The name of the file to write to.
	 */
	public void simulationStateToCSV(String filename) {
		arrayToCSV(filename, crystalMap);
	}
	
	/**
	 * Prints stats for the crystals in the simulation in the form of a CSV file.
	 * 
	 * @param filename The name of the file to write to.
	 */
	public void printCrystalStats(String filename) {
		try {
			File out = new File(filename);
			BufferedWriter output = new BufferedWriter(new FileWriter(out));
			finalizeCrystals();
			
			output.write("id,edge,init_area,final_area,init_perim,final_perim,#sticking_ops,%stuck_arrived_in_zone\n");
			for (short cId : crystals.keySet()) {
				Crystal c = crystals.get(cId);
				output.write(c.id + "," + c.edge + "," + c.initArea + "," + c.finalArea + "," + c.initPerim + "," + c.finalPerim);
				output.write("," + c.numStickingOps + "," + (1.0 * c.numArrived)/(c.numArrived + c.numStuck)+ "\n");
			}
			output.close();
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}
	}
	
	/**
	 * Writes stats for the monomers that have been used in the simulation to a file.
	 * Prints the time in existence and distance traveled for every monomer.
	 * 
	 * @param filename The name of the file to write to.
	 */
	public void printMonomerStats(String filename) {
		try {
			File out = new File(filename);
			BufferedWriter output = new BufferedWriter(new FileWriter(out));
			finalizeCrystals();
			
			output.write("time,dist\n");
			for (Monomer m : monomers) {
				output.write(m.timeExist + "," + m.distT + "\n");
			}
			output.close();
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}
	}
	
	/**
	 * Writes the contour map to a file in the form of a CSV file. Every entry is the
	 * number of MC cycles during which a monomer was on a particular pixel.
	 * 
	 * @param filename The name of the file to write to.
	 */
	public void printContourMap(String filename) {
		arrayToCSV(filename, contourDiagram);
	}
	
	/////////////////////////////////////////////////////////////////////////////////////
	////                          INNER STATIC CLASSES                               ////
	/////////////////////////////////////////////////////////////////////////////////////
	
	// Represents a single crystal in the MC simulation
	private static class Crystal {
		short id;           // the crystal's unique id
		boolean edge;       // true if this crystal's capture zone goes off the edge of the image
		int initArea;       // the pre-simulation area of the crystal
		int finalArea;      // the post-simulation area of the crystal
		int initPerim;      // the pre-simulation perimeter of the crystal
		int finalPerim;     // the post-simulation area of the crystal
		int numStickingOps; // the number of cycles a monomer has been adjacent to this crystal
		int numArrived;     // the # of monomers that arrived in this crystal's capture zone
		int numStuck;       // of those that arrived in this capture zone, the # that actually stuck
		
		// Creates a new crystal with the given id
		public Crystal(short id) {
			this.id = id;
			edge = false;
			initArea = 0;
			finalArea = 0;
			initPerim = 0;
			numStickingOps = 0;
			numArrived = 0;
			numStuck = 0;
		}
		
		@Override
		public int hashCode() {
			return Short.hashCode(id);
		}
		@Override
		public boolean equals(Object o) {
			if (!(o instanceof Crystal))
				return false;
			return ((Crystal)o).id == this.id;
		}		
	}
	
	/**
	 * Represents a single monomer (random walker).
	 */
	private static class Monomer {
		double timeExist;  // time in existence
		int distT;  	   // distance traveled
		
		// Creates a new monomer
		public Monomer() {
			timeExist = 0.0;
			distT = 0;
		}
		
		// Updates the time for the monomer at the end of a cycle
		// using the equation:
		//
		// t = t + ln(1/u') / R_n
		//
		// where R_n is the sum of the probabilities of all possible
		// moves for the monomer at this cycle, and u' is a number
		// chosen uniformly from 0 to 1
		public void updateTime(double R_n) {
			timeExist += Math.log(1 / random.nextDouble()) / R_n;
		}
	}
	
	/////////////////////////////////////////////////////////////////////////////////////
	////                       STATIC HELPER I/O METHODS                             ////
	/////////////////////////////////////////////////////////////////////////////////////
	
	public static void arrayToCSV(String filename, short[][] pixels) {
		try {
			File out = new File(filename);
			BufferedWriter output = new BufferedWriter(new FileWriter(out));
			int rows = pixels.length;
			int columns = pixels[0].length;
			for (int j = 0; j < rows; j++) {
				for (int i = 0; i < columns; i++) {
					if (i == columns - 1) {
						// last pixel on the row; print a newline instead of ','
						output.write(pixels[j][i] + "\n");
					}else {
						output.write(pixels[j][i] + ",");
					}
				}
			}
			output.close();
		} catch (IOException e) {
			System.err.println("Caught IOException: " + e.getMessage());
			System.exit(1);
		}
	}
	
	public static short[][] arrayFromCSV(String filename) {
		short[][] crystalMap = null;
		try (BufferedReader r = new BufferedReader(new FileReader(filename))) {
			// Get the height in pixels of the image
			int rows = getLineCount(filename);
			String line;
			int lineNumber = 0;
		    while ((line = r.readLine()) != null) {
		       String[] rowVals = line.split(",");
		       if (crystalMap == null) {
		    	   // this is the first line we're reading; get the width and
		    	   // initialize the crystal map
		    	   int columns = rowVals.length;
		    	   crystalMap = new short[rows][columns];
		       }
		       
		       for (int i = 0; i < crystalMap[0].length; i++) {
		    	   crystalMap[lineNumber][i] = Short.parseShort(rowVals[i]); 
		       }
		       lineNumber++;
		    }
		} catch (IOException e){
			 System.err.println("Caught IOException: " + e.getMessage());
			 System.exit(1);
		}
		return crystalMap;
	}
	
	public static short[][] computeTesselation(short[][] crystalMap) {
		short[][] tesselation = new short[crystalMap.length][crystalMap[0].length];
		for (int i = 0; i < crystalMap.length; i++) {
			for (int j = 0; j < crystalMap[0].length; j++) {
				short closestCrystal = computeClosestCrystal(crystalMap, i, j);
				tesselation[i][j] = closestCrystal;
			}
			System.out.println(i + "/" + crystalMap.length);
		}
		return tesselation;
	}
	
	private static short computeClosestCrystal(short[][] crystalMap, int row, int column) {
		if (crystalMap[row][column] != 0)
			// we are on a crystal (so closest crystal is this one)
			return crystalMap[row][column];
		
		short closestCrystal = -1;
		double minDist = Double.MAX_VALUE;
		for (int i = 0; i < crystalMap.length; i++) {
			for (int j = 0; j < crystalMap[0].length; j++) {
				if (crystalMap[i][j] == 0) 
					continue;  // this point is not a crystal
				double xDist = row - i;
				double yDist = column - j;
				// we can skip the sqrt sine we're just interested in relative dists
				double dist = xDist * xDist + yDist * yDist;
				if (dist < minDist) {
					minDist = dist;
				    closestCrystal = crystalMap[i][j];
				}
			}
		}
		assert (closestCrystal != -1) : "Closest crystal not found";
		
		return closestCrystal;
	}
	
	/**
	 * Returns the line count for the given file.
	 */
	public static int getLineCount(String filename) throws IOException {
		LineNumberReader  l = new LineNumberReader(new FileReader(new File(filename)));
		l.skip(Long.MAX_VALUE);
		int numLines = l.getLineNumber();
		l.close();
		return numLines + 1;	
	}
	
	/////////////////////////////////////////////////////////////////////////////////////
	////                                   MAIN                                      ////
	/////////////////////////////////////////////////////////////////////////////////////
	
	public static void main(String[] args) {
		MCSimulator walker = new MCSimulator("src/results.csv");
		//arrayToCSV("tessellation.csv", walker.tessellation);
		walker.simulationStateToCSV("initialState.csv");
		for (int i = 0; i < 30000; i++) {
			/*if (i % 250 == 0) {
				walker.simulationStateToCSV("output );
			}*/
			walker.doRandomWalk();
			System.out.println("walk " + (i + 1));
		}
		walker.simulationStateToCSV("finalState.csv");
		walker.printCrystalStats("finalStats.csv");
		walker.printMonomerStats("monomerStats.csv");
		walker.printContourMap("contourMap.csv");
	}
}
