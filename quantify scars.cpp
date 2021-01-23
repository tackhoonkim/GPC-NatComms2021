#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
//This version can find exact match only
//For analyzing paired end sequencing from MGH

using namespace std;

int match(string read, string dict[])
{
	//read: NGS read, dict: list of index/ampliconseq
	//size: number of different DNA seq to look up
	//length: length of DNA seq to look up
	//key_location: where index/ampliconseq should be found
	int size=0;
	while(dict[size]!="end")
		size++;

	for (int i=0;i<size;i++)
	{
		if (dict[i]==read.substr(20,6))
		{
			return i;
		}
	}
	return -1;	//no match in dictionary
}


int main() {
	int samples=18;		//number of samples in NGS analysis
	int loci=4;			//number of genomic loci in NGS analysis

/*	string index[samples+1]={"ACGTTG","AGACAG","AGAGCT",
					"CCAGTA","CTATCG","GAAGCA",
					"GAGGAT","GGATAC","GGCATA",
					"GTCTCA","GTGAGA","TAAGGC",
					"TGTGCA","TTCCGT","TTGCCA",
					"CACCAA","GGTCTT","TCTCCT","end"};
*/	string ampliconseq[loci+1]={"ATAACT","CCGAAG","CCCCCA","GGAGTA","end"};	//LoxP, FRT, AttP, AttL

	int counts[samples*loci];
	for (int i=0; i<samples*loci;i++)
		counts[i]=0;

	string read_sequence;	//parsing each line in seq file

	//Reads seq file
	ifstream seq[samples];
	ofstream results;
	results.open("results.txt");
	results<<"\tstage1\tstage2\tstage3\tstage4\n";
	int amplicon_found;
	
	for (int i=0; i<samples;i++)
	{
		seq[i].open("scar-"+to_string(i)+"-F.fastq");
		int counts[loci];
		for (int j=0; j<loci;j++)
			counts[j]=0;
		int counter=0;
		while (getline(seq[i],read_sequence))
		{
			counter++;
			if(counter%4==2)
			{
				read_sequence.erase(read_sequence.find_last_not_of(" \r\n\t")+1);		
				amplicon_found=match(read_sequence,ampliconseq);	//identify index

				if (amplicon_found>-1)
					counts[amplicon_found]++;
			}
		}
		seq[i].close();
		seq[i].open("scar-"+to_string(i)+"-R.fastq");
		counter=0;
		while (getline(seq[i],read_sequence))
		{
			counter++;
			if(counter%4==2)
			{
				read_sequence.erase(read_sequence.find_last_not_of(" \r\n\t")+1);		
				amplicon_found=match(read_sequence,ampliconseq);	//identify index

				if (amplicon_found>-1)
					counts[amplicon_found]++;
			}
		}
		results<<i<<"\t";
		for (int k=0; k<loci;k++)
			results<<counts[k]<<"\t";
		results<<"\n";
		seq[i].close();	
	}
	results.close();
	return 0;
}
