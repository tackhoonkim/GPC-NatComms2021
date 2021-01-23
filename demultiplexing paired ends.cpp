#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
//This version can find exact match only
//For analyzing paired end sequencing from MGH

using namespace std;

int match(string read, string dict[],int size, int key_location)
{
	//read: NGS read, dict: list of index/ampliconseq
	//key_location: where index/ampliconseq should be found
	for (int i=0; i<size;i++)
	{
		if (read.find(dict[i])==key_location)
			return i;
	}
	return -1;	//no match in dictionary
}

int main() {
	int samples=18;		//number of samples in NGS analysis
	int loci=5;			//number of genomic loci in NGS analysis

	string index[samples]={"ACGTTG","AGACAG","AGAGCT",
					"CCAGTA","CTATCG","GAAGCA",
					"GAGGAT","GGATAC","GGCATA",
					"GTCTCA","GTGAGA","TAAGGC",
					"TGTGCA","TTCCGT","TTGCCA",
					"CACCAA","GGTCTT","TCTCCT"};
	string ampliconseq[loci]={"GGTTTT","AAAAAG","GTTGGG","TTTCAC","CCTGCA"};	//APC, MLH1, SMAD4, TP53 PCR amplicon seq adjacent to index-less primer.
	string names[loci]={"APC","MLH1","SMAD4","TP53","scar"};
	
	string read_sequence;	//parsing each line in fastq file
		
	ofstream read1[samples*loci];
	for (int i=0;i<loci;i++)
	{
		for (int j=0; j<samples;j++)
			read1[i*samples+j].open(names[i]+"-"+to_string(j)+"-F.fastq");
	}
	ofstream read2[samples*loci];
	for (int i=0;i<loci;i++)
	{
		for (int j=0; j<samples;j++)
			read2[i*samples+j].open(names[i]+"-"+to_string(j)+"-R.fastq");
	}

	//Reads Fastq file
	ifstream fastqf, fastqr;
	fastqf.open("NGS-F.fastq");
	fastqr.open("NGS-R.fastq");
	string entry1, entry2;	//entire FASTQ entry
	string seq1, seq2;		//Sequence in each FASTQ entry
	
	while (getline(fastqf,read_sequence))
	{
		int index_found, amplicon_found;
		entry1=read_sequence+"\n";
		entry2="";
		for (int i=0; i<3; i++)
		{
			getline(fastqf,read_sequence);
			entry1=entry1+read_sequence+"\n";
			if (i==0)
				seq1=read_sequence;
		}
		for (int i=0; i<4; i++)
		{
			getline(fastqr,read_sequence);
			entry2=entry2+read_sequence+"\n";
			if (i==1)
				seq2=read_sequence;
		}		
		
		seq1.erase(seq1.find_last_not_of(" \r\n\t")+1);
		seq2.erase(seq2.find_last_not_of(" \r\n\t")+1);
		index_found=match(seq1,index,sizeof(index)/sizeof(index[0]),0);	//identify index
		if (index_found>-1)
		{
			int found1, found3;
			found1=entry1.find("\n");	//position of first new line
			found3=entry1.find("\n",found1+2);	//second
			found3=entry1.find("\n",found3+2);	//third
			entry1=entry1.substr(0,found1+1)+
					entry1.substr(found1+7,found3-found1-6)+
					entry1.substr(found3+7);		//trim multiplex adapter for sequence and Phred score
			amplicon_found=match(seq2,ampliconseq,sizeof(ampliconseq)/sizeof(ampliconseq[0]),20);
		}
		else
		{
			index_found=match(seq2,index,sizeof(index)/sizeof(index[0]),0);
			int found1, found3;
			found1=entry2.find("\n");	//position of first new line
			found3=entry2.find("\n",found1+2);	//second
			found3=entry2.find("\n",found3+2);	//third
			entry2=entry2.substr(0,found1+1)+
					entry2.substr(found1+7,found3-found1-6)+
					entry2.substr(found3+7);		//trim multiplex adapter for sequence and Phred score
			amplicon_found=match(seq1,ampliconseq,sizeof(ampliconseq)/sizeof(ampliconseq[0]),20);
		}
		//cout<<amplicon_found<<"\n";
		if (index_found>-1 && amplicon_found>-1)
		{
			read1[amplicon_found*samples+index_found]<<entry1;
			read2[amplicon_found*samples+index_found]<<entry2;
		}
	}
	for (int i=0;i<loci*samples;i++)
	{
		read1[i].close();
		read2[i].close();
	}
	return 0;
}
