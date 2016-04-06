import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;
import java.util.Map.Entry;

public class GlobalDDS {
	List<List<Double>> table;
	List<String> decision;
	List<String> headers;
	int rowCount;
	int attributeCount;
	Map<Integer, List<Double>> cutpoints;
	Map<Integer, String> decisionAssociation;
	List<List<Integer>> decisionSets;
	Map<Integer, Map<Double, Integer>> uniqueFreqCount;
	static Scanner scanner;
	int[] attrK;
	GlobalDDS() {
		table = new ArrayList<List<Double>>();
		decision = new ArrayList<String>();
		cutpoints = new HashMap<Integer, List<Double>>();
		decisionAssociation = new HashMap<Integer, String>();
		decisionSets = new ArrayList<List<Integer>>();
	}

	public static double round(double value, int places) {
		if (places < 0) throw new IllegalArgumentException();

		BigDecimal bd = new BigDecimal(value);
		bd = bd.setScale(places, RoundingMode.HALF_UP);
		return bd.doubleValue();
	}

	static double log(double x, int base)
	{
		return (Math.log(x) / Math.log(base));
	}

	public void readFile() {
		boolean fileFound = false;
		String input = "";
		do {
			System.out.println("Please enter the input filename containing the dataset to be discretized.\n Enter 4 to exit");
			input = scanner.nextLine();
			try {
				BufferedReader br = new BufferedReader(new FileReader(input));
				int currentInput = br.read();
				StringBuilder sb = new StringBuilder();
				List<Double> row = new ArrayList<Double>();
				boolean headerRow = false;
				int currAttrCnt = 0;
				while (currentInput != -1) {
					char currentChar = (char) currentInput;
					if (currentChar == '!') {
						while (currentInput != -1 && currentChar != '\n') {
							currentInput = br.read();
							currentChar = (char)currentInput;
						}
						currentInput = br.read();
						continue;
					}
					if (currentChar == '<') {
						while (currentInput != -1 && currentChar != '>') {
							currentInput = br.read();
							currentChar = (char)currentInput;
						}
						currentInput = br.read();
						continue;
					}
					
					if (currentChar == '+' ||
						currentChar == '-' ||
						currentChar == '.' ||
						currentInput >= 48 && currentInput <= 57 || //numbers
						currentInput >= 65 && currentInput <= 90 || //lowercase chars
						currentInput >= 97 && currentInput <= 122) { // uppercase chars
						sb.append(currentChar);
						currentInput = br.read();
						continue;
					}
					
					if (currentChar == '[') {
						headerRow = true;
						headers = new ArrayList<String>();
						currentInput = br.read();
						continue;
					}
					
					if (currentChar == ']') {
						headerRow = false;
						attributeCount = headers.size() - 1;
						currentInput = br.read();
						continue;
					}
					
					if (currentChar == ' ' ||
						currentChar == '\t' ||
						currentChar == '\n') {
						if (sb.length() == 0) {
							currentInput = br.read();
							continue;
						}
						if (headerRow) {
							headers.add(sb.toString());
						} else {
							if (currAttrCnt == attributeCount) {
								table.add(new ArrayList<Double>(row));
								row.clear();
								decision.add(sb.toString());
								sb.setLength(0);
								currAttrCnt = 0;
							} else {
								double d = Double.parseDouble(sb.toString());
								row.add(d);
								currAttrCnt++;
							}
						}
						
						while (currentChar == ' ') {
							currentInput = br.read();
							currentChar = (char) currentInput;
						}
						sb.setLength(0);
						continue;
					}
					currentInput = br.read();
				}
				br.close();
				fileFound = true;
				rowCount = table.size();
				attributeCount = table.get(0).size();
				Map<String, List<Integer>> temp = new HashMap<String, List<Integer>>();
				for (int i = 0; i < rowCount; i++) {
					decisionAssociation.put(i, decision.get(i));
					if (temp.containsKey(decision.get(i))) {
						temp.get(decision.get(i)).add(i);
					} else {
						List<Integer> list = new ArrayList<Integer>();
						list.add(i);
						temp.put(decision.get(i), list);
					}
				}

				for (Entry<String, List<Integer>> e : temp.entrySet()) {
					decisionSets.add(e.getValue());
				}
			} catch (FileNotFoundException fe) {
				System.out.println("Input file "+ input +" not found!!!" + fe);
			} catch (Exception e) {
				System.out.println(e);
			}
		} while(!input.equals("4") && !fileFound);
	}

	public void solveByEqualFrequency() {
		attrK = new int[attributeCount];
		Arrays.fill(attrK, 2);
		List<List<String>> cutpointTable = constructInitialTableEF();
		while (!isDataConsistent(cutpointTable)) {
			cutpointTable = constructCutpointTableEF(cutpointTable);
		}
		mergeAndPrint();
	}
	
	public void solveByEqualInterval() {
		attrK = new int[attributeCount];
		Arrays.fill(attrK, 2);
		List<List<String>> cutpointTable = constructInitialTableEI();
		while (!isDataConsistent(cutpointTable)) {
			cutpointTable = constructCutpointTableEI(cutpointTable);
		}
		mergeAndPrint();
	}

	public void solveConditionalEntropy() {
		List<List<String>> cutpointTable = constructInitialTableCE();
		while (!isDataConsistent(cutpointTable)) {
			cutpointTable = constructCutpointTableCE(cutpointTable);
		}
		mergeAndPrint();
	}

	private void mergeAndPrint() {
		try {
		PrintWriter cutpointFile = new PrintWriter("test.int", "UTF-8");
		cutpointFile.println("Cutpoints BEFORE merging\n\n");
		for (int j = 0; j < attributeCount; j++) {
			cutpointFile.println("Cutpoints for attribute " + headers.get(j) + " are:");
			List<Double> attrCutpoints = cutpoints.get(j);
			for (double cutpoint : attrCutpoints)
				cutpointFile.print(cutpoint + " ");
			cutpointFile.println();
		}
		cutpointFile.println("\n\n");
		
		//merging
		Map<Integer, List<Double>> finalCutpoints = new HashMap<Integer, List<Double>>(cutpoints);
		//try removing cutpoints for each attribute
		for (int j = 0; j < attributeCount; j++) {
			List<Double> attrCutpoints = cutpoints.get(j);
			for (int i = 0; i < attrCutpoints.size(); i++) {
				Double cutpoint = attrCutpoints.get(i);
				Map<Integer, List<Double>> currCutpoints = new HashMap<Integer, List<Double>>(cutpoints);
				for (Entry<Integer, List<Double>> e : cutpoints.entrySet()) {
					currCutpoints.put(e.getKey(), new ArrayList<Double>(e.getValue()));
				}
				currCutpoints.get(j).remove(cutpoint);
				List<List<String>> modifiedCutpointTable = buildCutpointTable(currCutpoints);
				if (isDataConsistent(modifiedCutpointTable)) {
					finalCutpoints.get(j).remove(cutpoint);
				}
			}
		}
		
		
			
			PrintWriter dataFile = new PrintWriter("test.data", "UTF-8");
			
			cutpointFile.println("Cutpoints AFTER merging\n\n");
			//print cutpoints
			for (int j = 0; j < attributeCount; j++) {
				cutpointFile.println("Cutpoints for attribute " + headers.get(j) + " are:");
				List<Double> attrCutpoints = finalCutpoints.get(j);
				for (double cutpoint : attrCutpoints)
					cutpointFile.print(cutpoint + " ");
				cutpointFile.println();
			}
			
			//print final table
			List<List<String>> finalTable = buildCutpointTable(finalCutpoints);
			
			dataFile.print("< ");
			for (int i = 0; i < attributeCount; i++) {
				dataFile.print("a ");
			}
			dataFile.print("d >");
			dataFile.println();
			
			dataFile.print("[ ");
			for (int i = 0; i <= attributeCount; i++) {
				dataFile.print(headers.get(i) + " ");
			}
			dataFile.print("]");
			dataFile.println();
			for (int i = 0; i < rowCount; i++) {
				for (int j = 0; j < attributeCount; j++) {
					dataFile.print(finalTable.get(i).get(j) + " ");
				}
				dataFile.print(decision.get(i));
				dataFile.println();
			}
			cutpointFile.close();
			dataFile.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private List<List<String>> constructCutpointTableEF(
			List<List<String>> cutpointTable) {
		int worstEntropyAttr = getWorstBlockEntropy(cutpointTable);
		attrK[worstEntropyAttr]++;
		calcEFCutpoints();
		List<List<String>> res = buildCutpointTable(cutpoints); 
		return res;
	}
	
	private List<List<String>> constructCutpointTableEI(
			List<List<String>> cutpointTable) {
		int worstEntropyAttr = getWorstBlockEntropy(cutpointTable);
		attrK[worstEntropyAttr]++;
		calcEICutpoints();
		List<List<String>> res = buildCutpointTable(cutpoints); 
		return res;
	}


	private List<List<String>> constructCutpointTableCE(
			List<List<String>> cutpointTable) {
		int worstEntropyAttr = getWorstBlockEntropy(cutpointTable);
		Map<String, List<Integer>> common = new HashMap<String, List<Integer>>();
		for (int i = 0; i < rowCount; i++) {
			if (common.containsKey(cutpointTable.get(i).get(worstEntropyAttr))) {
				common.get(cutpointTable.get(i).get(worstEntropyAttr)).add(i);
			} else {
				List<Integer> indexes = new ArrayList<Integer>();
				indexes.add(i);
				common.put(cutpointTable.get(i).get(worstEntropyAttr), indexes);
			}
		}

		double minEntropy = Double.MAX_VALUE;
		double minCutpoint = -1;
		//get another cutpoint from the subtables of the worst entropy attribute
		for (Entry<String, List<Integer>> e : common.entrySet()) {
			List<Double> subtable = new ArrayList<Double>();
			List<String> subDecision = new ArrayList<String>();
			List<Integer> entries = e.getValue();
			Set<Double> unique = new HashSet<Double>();
			double minVal = Double.MAX_VALUE;
			double maxVal = Double.MIN_VALUE;
			for (int i = 0; i < entries.size(); i++) {
				subtable.add(table.get(entries.get(i)).get(worstEntropyAttr));
				unique.add(table.get(entries.get(i)).get(worstEntropyAttr));
				subDecision.add(decision.get(entries.get(i)));
				minVal = Math.min(minVal, table.get(entries.get(i)).get(worstEntropyAttr));
				maxVal = Math.max(maxVal, table.get(entries.get(i)).get(worstEntropyAttr));
			}
			if (unique.size() == 0)
				continue;
			List<Double> uniqueList = new ArrayList<Double>(unique);
			int len = uniqueList.size();
			for (int i = 1; i < len; i++) {
				double currCutpoint = (uniqueList.get(i - 1) + uniqueList.get(i)) / 2.0;
				currCutpoint = round(currCutpoint, 2);
				List<String> modifiedAttr = new ArrayList<String>();
				for (int j = 0; j < subtable.size(); j++) {
					if (subtable.get(j) < currCutpoint) {
						modifiedAttr.add(minVal + ".." + currCutpoint);
					} else {
						modifiedAttr.add(currCutpoint + ".." + maxVal);
					}
				}
				double currEntropy = getAttributeEntropy(modifiedAttr, subDecision);
				if (currEntropy < minEntropy) {
					minEntropy = currEntropy;
					minCutpoint = currCutpoint;
				}
			}
		}

		//add the new cutpoint to the original list of cutpoints
		minCutpoint = round(minCutpoint, 2);
		cutpoints.get(worstEntropyAttr).add(minCutpoint);
		//create the new cutpoint table
		List<List<String>> result = buildCutpointTable(cutpoints);
		return result;
	}

	private int getWorstBlockEntropy(List<List<String>> cutpointTable) {
		//get the worst block entropy column
		double worstEntropy = Double.MIN_VALUE;
		int worstEntropyAttr = -1;
		for (int j = 0; j < attributeCount; j++) {
			List<String> column = new ArrayList<String>();
			for (int i = 0; i < rowCount; i++) {
				column.add(cutpointTable.get(i).get(j));
			}
			int cutpointCnt = cutpoints.get(j).size();
			double conditionalEntropy = getAttributeEntropy(column, decision);
			double blockEntropy = conditionalEntropy / (cutpointCnt + 1);
			if (blockEntropy > worstEntropy) {
				worstEntropy = blockEntropy;
				worstEntropyAttr = j;
			}
		}
		return worstEntropyAttr;
	}

	private List<List<String>> buildCutpointTable(Map<Integer, List<Double>> cutpoints) {
		List<List<String>> res = new ArrayList<List<String>>();
		for (int i = 0; i < rowCount; i++) {
			res.add(new ArrayList<String>());
		}

		for (int j = 0; j < cutpoints.size(); j++) {
			List<Double> cutpointList = new ArrayList<Double>(cutpoints.get(j));
			double minVal = Double.MAX_VALUE;
			double maxVal = Double.MIN_VALUE;
			for (int i = 0; i < rowCount; i++) {
				minVal = Math.min(minVal, table.get(i).get(j));
				maxVal = Math.max(maxVal, table.get(i).get(j));
			}
			cutpointList.add(0, minVal);
			cutpointList.add(maxVal);
			Collections.sort(cutpointList);
			for (int i = 0; i < rowCount; i++) {
				for (int cnt = 1; cnt < cutpointList.size(); cnt++) {
					double tableVal = table.get(i).get(j);
					if (tableVal >= cutpointList.get(cnt - 1) &&
							tableVal <= cutpointList.get(cnt)) {
						res.get(i).add(cutpointList.get(cnt - 1) + ".." + cutpointList.get(cnt));
						break;
					}
				}
			}
		}
		return res;
	}

	private void calcMinimumEntropyCE(Map<Integer, Set<Double>> cutpointTable) {
		//get the minimum entropy for each attribute by going through each cutpoint
		for (int j = 0; j < attributeCount; j++) {
			List<Double> attrCutpoints = new ArrayList<Double>(cutpointTable.get(j));
			double minimumEntropy = Double.MAX_VALUE;
			double minimumCutpoint = Double.MAX_VALUE;
			//no cutpoints for the attribute, implies all values are equal. Take the first one
			if (attrCutpoints.size() == 0) {
				cutpoints.get(j).add(table.get(0).get(j));
				continue;
			}
			for (int i = 0; i < attrCutpoints.size(); i++) {
				double cutpoint = attrCutpoints.get(i);
				double minVal = Double.MAX_VALUE;
				double maxVal = Double.MIN_VALUE;
				for (int row = 0; row < rowCount; row++) {
					minVal = Math.min(minVal, table.get(row).get(j));
					maxVal = Math.max(maxVal, table.get(row).get(j));
				}
				List<String> modifiedAttr = new ArrayList<String>();
				for (int row = 0; row < rowCount; row++) {
					if (table.get(row).get(j) < cutpoint) {
						modifiedAttr.add(minVal + ".." + cutpoint);
					} else {
						modifiedAttr.add(cutpoint + ".." + maxVal);
					}
				}
				double currEntropy = getAttributeEntropy(modifiedAttr, decision);
				if (currEntropy < minimumEntropy) {
					minimumEntropy = currEntropy;
					minimumCutpoint = cutpoint;
				}
			}
			cutpoints.get(j).add(minimumCutpoint);
		}
	}

	private List<List<String>> constructInitialTableEF() {
		//create initial set of cutpoints
		calcEFCutpoints();
		List<List<String>> res = buildCutpointTable(cutpoints); 
		return res;
	}
	
	private List<List<String>> constructInitialTableEI() {
		//create initial set of cutpoints
		calcEICutpoints();
		List<List<String>> res = buildCutpointTable(cutpoints); 
		return res;
	}

	private List<List<String>> constructInitialTableCE() {
		//create initial set of cutpoints
		Map<Integer, Set<Double>> cutpointTable = getCECutpointTable(); 
		calcMinimumEntropyCE(cutpointTable);
		List<List<String>> res = buildCutpointTable(cutpoints); 
		return res;
	}

	private void calcEFCutpoints() {
		cutpoints.clear();
		for (int j = 0; j < attributeCount; j++) {
			cutpoints.put(j, new ArrayList<Double>());
			int K = attrK[j];
			double minVal = Double.MAX_VALUE;
			double maxVal = Double.MIN_VALUE;
			double[] values = new double[rowCount];
			for (int i = 0; i < rowCount; i++) {
				minVal = Math.min(minVal, table.get(i).get(j));
				maxVal = Math.max(maxVal, table.get(i).get(j));
				values[i] = table.get(i).get(j);
			}
			Arrays.sort(values);
			int cutpointCount = K - 1;
			for (int i = 1; i <= cutpointCount; i++) {
				double cutpoint = (values[(rowCount / K) * i] + values[(rowCount / K) * i + 1]) / 2.0;
				cutpoint = round(cutpoint, 2);
				cutpoints.get(j).add(cutpoint);
			}
		}
	}

	private void calcEICutpoints() {
		cutpoints.clear();
		for (int j = 0; j < attributeCount; j++) {
			cutpoints.put(j, new ArrayList<Double>());
			int K = attrK[j];
			double minVal = Double.MAX_VALUE;
			double maxVal = Double.MIN_VALUE;
			for (int i = 0; i < rowCount; i++) {
				minVal = Math.min(minVal, table.get(i).get(j));
				maxVal = Math.max(maxVal, table.get(i).get(j));
			}
			double interval = (maxVal - minVal) / K;
			for (int i = 1; i < K; i++) {
				double currCutpoint = minVal + interval * i;
				currCutpoint = round(currCutpoint, 2);
				cutpoints.get(j).add(currCutpoint);
			}
		}
	}

	private Map<Integer, Set<Double>> getCECutpointTable() {
		Map<Integer, Set<Double>> cutpointTable = new HashMap<Integer, Set<Double>>();
		for (int j = 0; j < attributeCount; j++) {
			Set<Double> uniqueValues = new HashSet<Double>();
			for (int i = 0; i < rowCount; i++) {
				uniqueValues.add(table.get(i).get(j));
			}
			cutpointTable.put(j, new HashSet<Double>());
			cutpoints.put(j, new ArrayList<Double>());
			List<Double> uniqueList = new ArrayList<Double>(uniqueValues);
			Collections.sort(uniqueList);
			int len = uniqueList.size();
			for (int i = 1; i < len; i++) {
				cutpointTable.get(j).add((uniqueList.get(i - 1) + uniqueList.get(i)) / 2.0);
			}
		}
		return cutpointTable;
	}

	private boolean isDataConsistent(List<List<String>> cutpointTable) {
		List<List<Integer>> aStar = new ArrayList<List<Integer>>();
		Map<String, List<Integer>> map = new HashMap<String, List<Integer>>();
		for (int i = 0; i < rowCount; i++) {
			StringBuilder sb = new StringBuilder();
			for (int j = 0; j < attributeCount; j++) {
				sb.append(cutpointTable.get(i).get(j));
			}
			String str = sb.toString();
			if (map.containsKey(str)) {
				map.get(str).add(i);
			} else {
				List<Integer> list = new ArrayList<Integer>();
				list.add(i);
				map.put(str, list);
			}
		}
		for (Entry<String, List<Integer>> e : map.entrySet()) {
			aStar.add(e.getValue());
		}
		for (int i = 0; i < decisionSets.size(); i++) {
			List<Integer> set = decisionSets.get(i);
			for (int j = 0; j < set.size(); j++) {
				int dIndex = set.get(j);
				String decisionStr = decisionAssociation.get(dIndex);
				//go through each AStar set to find decision index
				for (int k = 0; k < aStar.size(); k++) {
					List<Integer> aStarSet = aStar.get(k);
					boolean found = false;
					for (int l = 0; l < aStarSet.size(); l++) {
						if (aStarSet.get(l) == dIndex) {
							found = true;
							break;
						}
					}
					if (found) {
						//go over again to check consistency
						for (int l = 0; l < aStarSet.size(); l++) {
							String currDecisionStr = decisionAssociation.get(aStarSet.get(l)); 
							if (!currDecisionStr.equals(decisionStr))
								return false;
						}
					}
				}
			}
		}
		return true;
	}

	private double getAttributeEntropy(List<String> attrColumn,
			List<String> decisionColumn) {
		class AttributeData {
			Map<String, Integer> associations;
			public int count;
			AttributeData() {
				associations = new HashMap<String, Integer>();
				count = 0;
			}
			public void addKey(String key) {
				if (associations.containsKey(key)) {
					int k = associations.get(key);
					associations.put(key, k + 1);
				} else {
					associations.put(key, 1);
				}
				count++;
			}			
		}

		Map<String, AttributeData> map = new HashMap<String, AttributeData>();
		int rowCnt = attrColumn.size();
		for (int i = 0; i < rowCnt; i++) {
			String key = attrColumn.get(i);
			if (map.containsKey(key)) {
				map.get(key).addKey(decisionColumn.get(i));
			} else {
				AttributeData attr = new AttributeData();
				attr.addKey(decisionColumn.get(i));
				map.put(key, attr);
			}
		}
		double res = 0.0;
		for (Entry<String, AttributeData> e : map.entrySet()) {
			Map<String, Integer> associations = e.getValue().associations;
			for (Entry<String, Integer> association : associations.entrySet()) {
				double temp = (double)association.getValue() / e.getValue().count; 
				res += (double)e.getValue().count / rowCnt * temp * log(temp, 2); 
			}	
		}
		return res * -1.0;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		 scanner = new Scanner(System.in);
		GlobalDDS discreteDataSets = new GlobalDDS();
		//method to read the file in to a data structure.
		discreteDataSets.readFile();

		System.out.println("Please select the discretization algorithm from the following three options\n 1. Equal Interval Width\n 2. "
				+ "Equal frequency per interval\n 3. Conditional Entropy \n OR \n Enter 4 to exit \n Enter 1 , 2 , 3 or 4 below");
		boolean algorithmSelected = false;
		while(!algorithmSelected){
			int optionSelected = scanner.nextInt();
			switch(optionSelected){
			case 1:
				// call function for equal interval width
				algorithmSelected = true;
				discreteDataSets.solveByEqualInterval();
				break;

			case 2:
				// call function for equal frequency per interval
				algorithmSelected = true;
				discreteDataSets.solveByEqualFrequency();
				break;
			case 3:
				// call function conditional entropy
				algorithmSelected = true;
				discreteDataSets.solveConditionalEntropy();
				break;
			case 4:
				return;

			default:
				System.out.println("Invalid option");
				System.out.println("Please select the discretization algorithm from the following three options\n 1. Equal Interval Width\n 2. "
						+ "Equal frequency per interval\n 3. Conditional Entropy \n OR \n Enter 4 to exit \n Enter 1 , 2 , 3 or 4 below");
			}
		}
		scanner.close();
	}
}
