package csNlrManuscript;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

public class IntegratedDomains {

	
	
	
	/**
	 * 
	 * Print out a list of NLR proteins together with NB-ARC and integrated domains. LRRs are excluded as well as AAA and NACHT. The latter two always come up when you have an NB-ARC domain.
	 * 
	 * 
	 * @param evalue_threshold
	 * 			e-value threshold for parsing the hmmer file. For the manuscript, we used 1E-5.
	 * @param inputNLRList
	 * 			The inputNLRList is the list of NLR genes. One identifyer per line. "NLR_genes.txt"
	 * @param inputHMMER
	 * 			The inputHMMER is the output of HMMER run on all RefSeq v1.0 genes. "iwgsc_refseq1.0_All_REPR_PROTEIN.hmmer.txt"
	 * @throws IOException
	 * 			input files not present or lack of access rights.
	 */
	public static void generateIDList(double evalue_threshold, File inputNLRList, File inputHMMER )throws IOException{
		Hashtable<String, Vector<String>> nlrs = new Hashtable<String, Vector<String>>();
		BufferedReader in = new BufferedReader(new FileReader(inputNLRList));
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			
			nlrs.put(inputline, new Vector<String>());
		}
		in.close();

		in = new BufferedReader(new FileReader(inputHMMER));
		String currentID= "";
		int[] covered = new int[0];
		Vector<String> domains = new Vector<String>();
		boolean nbarc_only = true;
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			if(inputline.startsWith("#")){
				continue;
			}
			String[] split = inputline.split("\\s+");
			String id = split[3].split("\\.")[0];
			if(nlrs.containsKey(id)){
			
				String domain = split[0];
				if(domain.startsWith("AAA")){
					continue;
				}
				if(domain.startsWith("LRR")){
					continue;
				}
				if(domain.startsWith("NACHT")){
					continue;
				}
				double evalue = Double.parseDouble(split[6]);
				int start = Integer.parseInt(split[19]);
				int end = Integer.parseInt(split[20]);
				
				if( !id.equalsIgnoreCase(currentID)){
					
					if(domains.size() >0 && !nbarc_only){
					
						System.out.print(currentID + "\t");
						
						Collections.sort(domains, new Comparator<String>(){public int compare(String s1, String s2){
							String[] split1 = s1.split("\t");
							String[] split2 = s2.split("\t");
							int i1 = Integer.parseInt(split1[1]);
							int i2 = Integer.parseInt(split2[1]);
							
							if( i1<i2 ){return -1;}
							if( i1>i2 ){return 1;}
							return 0;}});
						
						String s = "";
						for(Enumeration<String> myenum = domains.elements(); myenum.hasMoreElements();){
							String[] splitDomain = myenum.nextElement().split("\t");
							s = s + ";"+splitDomain[0]+":"+splitDomain[1] + "-"+splitDomain[2];
						}
						
						System.out.println(s.substring(1));
					}
					int l = Integer.parseInt(split[5]);
					currentID = id;
					covered = new int[l];
					domains = new Vector<String>();
					nbarc_only = true;
				}
				
				
				
				if(evalue < evalue_threshold){
					boolean add = true;
					for( int i = start-1; i< end; i++ ){
						if( covered[i] >0){
							add=false;
							break;
						}
					}
					if(add){
						domains.add(domain + "\t" + start + "\t" + end);
						for( int i = start-1; i< end; i++){
							covered[i]++;
						}
						if(!domain.equalsIgnoreCase("NB-ARC") && !domain.startsWith("LRR")){
							nbarc_only = false;
						}
					}
					
					
					
					
				}
				
				
			}
			
		}

		in.close();
		
	}
	
	}
