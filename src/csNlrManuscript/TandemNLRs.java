package csNlrManuscript;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Vector;

public class TandemNLRs {

	
	
		/**
		 * 
		 * print out pairs of genes that are in a "head-to-head" arrangement and a distance of less than "maxDistance".
		 * 
		 * The printout text can be added to the iTOL templade https://itol.embl.de/help/dataset_connections_template.txt
		 * 
		 * connectors are printed in gray except for those having one NLR in the connector clade F.
		 * 
		 * 
		 * @param maxDistance
		 * 			The maximum distance between NLR pairs. For the manuscript, we used 50,000.
		 * @param inputConnectorClade
		 * 			A list of NLRs that are in clade F
		 * @param inputNLRList
		 * 			A list of NLR genes
		 * @param highConfGff
		 * 			GFF of high confidence genes of RefSeq annotation v1.0
		 * @param lowConfGff
		 * 			GFF of low confidence genes of RefSeq annotation v1.0
		 * @throws IOException
		 */
		public static void getTandems( int maxDistance,File inputConnectorClade,File inputNLRList,File highConfGff, File lowConfGff)throws IOException{
			
			
			HashSet<String>connectors = new HashSet<String>();
			BufferedReader in = new BufferedReader(new FileReader(inputConnectorClade));
			for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
				connectors.add(inputline);
			}

			in.close();
			
			
			
			/*
			 * select NLR genes
			 */
			
			
			HashSet<String> usable = new HashSet<String>();
			 in = new BufferedReader(new FileReader(inputNLRList));

			for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
				usable.add(inputline);
			}

			in.close();

			
			/*
			 * Load HC NLRs
			 * 
			 */
			Vector<String>v = new Vector<String>();
									
			in = new BufferedReader(new FileReader(highConfGff));
			for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
				//chr1A	IWGSC_March2017	mRNA	7185497	7186507	60	+	.	ID=TraesCS1A01G012600.1;Parent=TraesCS1A01G012600;primconf=HC;secconf=HC1

				String[] split = inputline.split("\t");
				if(!split[2].equalsIgnoreCase("gene")){
					continue;
				}
				String id = inputline.split("ID=")[1].split(";")[0];
				
				if(!usable.contains(id)){
					continue;
				}
				
				v.add(inputline);
				
				
			}
			in.close();
			
			/*
			 * Load LC NLRs
			 * 
			 */
			in = new BufferedReader(new FileReader(lowConfGff));
			for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
				//chr1A	IWGSC_March2017	mRNA	7185497	7186507	60	+	.	ID=TraesCS1A01G012600.1;Parent=TraesCS1A01G012600;primconf=HC;secconf=HC1

				String[] split = inputline.split("\t");
				if(!split[2].equalsIgnoreCase("gene")){
					continue;
				}
				String id = inputline.split("ID=")[1].split(";")[0];
				if(!usable.contains(id)){
					continue;
				}
				v.add(inputline);
			}
			in.close();

			
			Collections.sort(v, new Comparator<String>(){public int compare(String s1, String s2){
				String[] split1 = s1.split("\t");
				String[] split2 = s2.split("\t");
				
				int i1 = split1[0].compareTo(split2[0]);
				if( i1!=0 ){
					return i1;
				}
				
				int start1 = Integer.parseInt(split1[3]);
				int start2 = Integer.parseInt(split2[3]);
				if( start1<start2){
					return -1;
				}
				if( start2 < start1){
					return 1;
				}
				
				return 0;	
			}
			});
			
		
			
			
			for(int i = 0; i< v.size()-1; i++){
				String[] s1 = v.get(i).split("\t");
				int j = 1;
				String[] s2 = v.get(i+j).split("\t");
				int start1 = Integer.parseInt(s1[3]);
				int start2 = Integer.parseInt(s2[3]);
				
				while(j< v.size()-i && s1[0].equalsIgnoreCase(s2[0]) && start2-start1<maxDistance ){
					if(s1[0].equalsIgnoreCase(s2[0]) && s1[6].equalsIgnoreCase("-") && s2[6].equalsIgnoreCase("+")){
						String id1 = s1[8].split("ID=")[1].split(";")[0];
						String id2 = s2[8].split("ID=")[1].split(";")[0];
						
						String color = "#808080"; //a gray connector for standard
						if(connectors.contains(id1) || connectors.contains(id2)){
							color = "#2BCE48"; // a green connector for those coming out of the "connectorClade"
						}
						
						System.out.println(id1 + "," + id2 +",1,"+color+",");
					}
					j++;
					s2 = v.get(i+j).split("\t");
					start2 = Integer.parseInt(s2[3]);
				}
				
				
				
				
				
				
			}
			
			
		}
		

		
		
		
		
		
		
		
	
		
		
	
}
