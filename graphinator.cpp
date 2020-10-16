/*
GRAPHANATOR!
by Kim Ng

Some older code I found which I did at KU related to bayesembler. From 14/06/14. What it was for was creating a graph from a bam file to find all possible isoforms paths.

Currently looking at Splice Graph creation adding extra edges, occuring when edge starts from same segment from midway and goes to different locations
*/

/*---std library---*/
#include <iostream>
#include <sstream>
#include <algorithm>
#include <ctime> //for timing of application run
#include <stack>
#include <string>

/*---BamTools---*/
#include "api/BamReader.h"
#include "api/BamWriter.h"

/*---Boost library---*/
#include <boost/config.hpp>
#include <boost/unordered_map.hpp>//for EdgePoints
#include <boost/graph/adjacency_list.hpp> //graph is made from adjacency list
#include <boost/graph/graphviz.hpp> //used to display graph
#include <boost/graph/graph_utility.hpp> //for debugging ie print_graph
//#include <boost/log/trivial.hpp>
#include <boost/algorithm/string.hpp> //additional string manipulations used
#include <boost/bind.hpp>

//#include <boost/graph/filtered_graph.hpp> //used to filter graph on removed nodes
//#include <boost/graph/incremental_components.hpp> //
//#include <boost/pending/disjoint_sets.hpp> //
//#include <boost/progress.hpp>//messing around with this for progress bars to improve user experience

using namespace std;
using namespace BamTools;
using namespace boost;

time_t startTime, processTime, totalTime;

#define DEBUG//undocument for debug mode
class MyVertex
{
	public:
		//Core parts of the vertex's to be used later
		int startPosition;
		int endPosition;		
		//For naming in the graph output
		std::string name;		
};
class MyEdge
{
	public:		
};

/*---Custom Types---*/
typedef boost::unordered_map<std::pair<int, int>, int> EdgePointMap;
typedef boost::unordered_map<std::pair<int, int>, int> SegmentMap; //for GTF files
/*---Functions---*/

/*
Function: returns all of the files contents as a string
*/
std::string getFileContents(const char *filename)
{
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (in)
  {
    std::string contents;
    in.seekg(0, std::ios::end);
    contents.resize(in.tellg());
    in.seekg(0, std::ios::beg);
    in.read(&contents[0], contents.size());
    in.close();
    return(contents);
  }
}

/************************************************************************************************************
Function:
//Checks the conditions of the alignment read against the pre-inputed requirements. Currently checks if the alignment is on the + or - strand, that the reference ID matches the current chromosome and that the rea is only mapped to one location.
//to bypass check on strand use o and -1 for Chromosome Reference ID
************************************************************************************************************/
class RunConditions 
{
	public:
		char runStrand;
		string chromosomeReferenceName;
};
bool runConditionsCheck(BamAlignment alignment, char runStrand, int chromosomeReferenceID)
{
	char strand;
	uint32_t nhTag;//tag info requires specialized type		
	alignment.GetTag("NH",nhTag);//check Number of Hits 

	if(	(alignment.IsFirstMate() and alignment.IsReverseStrand()) or 
		(!alignment.IsFirstMate() and !alignment.IsReverseStrand())) 
	{
		strand = '+';
	} 
	else if((alignment.IsFirstMate() and !alignment.IsReverseStrand()) or 
			(!alignment.IsFirstMate() and alignment.IsReverseStrand()))  
	{
		strand = '-';
	}	
	if(nhTag<=1 and (strand==runStrand or runStrand=='o') and (alignment.RefID==chromosomeReferenceID or chromosomeReferenceID==-1)) 
	{
		return true;
	}
	else
	{//Otherwise skip the alignment
		return false;
	}
}
//Debug Function
int printCurrentRunTime(time_t &startTime, time_t &processTime, time_t &totalTime)
{
	#ifdef DEBUG
	double processRunTime = difftime(processTime,startTime);//previous time	
	time(&totalTime);
	double totalRunTime = difftime(totalTime,startTime);
	time(&processTime);
	std::cout << "Total runtime: "<< totalRunTime << "\tprocess runtime: " << totalRunTime-processRunTime << std::endl;
	#endif
	return 0;
}
int printCurrentRunTime(time_t &startTime, time_t &processTime, time_t &totalTime, std::string message)
{
	#ifdef DEBUG
	double processRunTime = difftime(processTime,startTime);//previous time	
	time(&totalTime);
	double totalRunTime = difftime(totalTime,startTime);
	time(&processTime);
	std::cout << "Total runtime: "<< totalRunTime << "\tprocess runtime: " << totalRunTime-processRunTime << "\t" << message << std::endl;
	#endif
	return 0;
}

int chromosomeReferenceIDFromReferenceName(RefVector references, std::string chromosomeReferenceName)
{//return id or -1 otherwise	
	for(int i=0; i<references.size();i++) 
	{	
		if(references[i].RefName==chromosomeReferenceName)
		{
			return i;
		}
	}
	return -1;
}

//Debug Function
class Segment {
	public:	
		Segment() 
		{
			this->position = std::make_pair(0,0);
		}				
		Segment(int startPosition, int length) 
		{
			this->position.first = startPosition; 
			this->position.second = startPosition + length - 1;
		}
		Segment(std::pair<int, int> position)
		{
			this->position = position;
		}
		inline const int& startPosition() const { return position.first; }
		inline const int& endPosition() const { return position.second; }
		inline const pair<int, int> positions() {return position; }
		int Length() { return (position.second - position.first +1 ); }			

		bool checkForOverlap(Segment comparisonSegment)
		{
			if(	this->startPosition() <= comparisonSegment.startPosition() and
				this->endPosition() >= comparisonSegment.startPosition() )
			{
				return true;
			}
			else if( comparisonSegment.startPosition() <= this->startPosition() and
					 comparisonSegment.endPosition() >= this->startPosition() )
			{
				return true;
			}
			else if( this->endPosition()+1 == comparisonSegment.startPosition() or
					 this->startPosition() == comparisonSegment.endPosition()+1 )
			{//just bordering each other also means its okay to merge
				return true;
			}
			else
			{
				return false;
			}
		}		
		bool extendEndPositionWithSegment(Segment segment)
		{	
			if(	segment.startPosition() <= this->endPosition() and 
				segment.startPosition() >= this->startPosition())
			{				
				this->position.second = max(this->endPosition(),segment.endPosition());
				return true;
			}
			else
			{
				return false;
			}
		}
		bool mergeIntoSegment(Segment segment)		
		{
			if(checkForOverlap(segment))
			{
				this->position.first = std::min(this->startPosition(), segment.startPosition());
				this->position.second = std::max(this->endPosition(), segment.endPosition());
				return true;
			}			
			else
			{
				return false;
			}
		}	
	private:
		std::pair <int, int> position;		
};
bool compareSegments(const Segment &lhs, const Segment &rhs)
{
	return lhs.startPosition() < rhs.startPosition();
}
class EdgePoint {
	public:	
		EdgePoint() {}	
		EdgePoint(int startPosition, int length) 
		{
			this->position.first = startPosition;
			this->position.second = startPosition + length + 1;
		}
		EdgePoint(std::pair<int, int> position, int count)
		{
			this->position.first = position.first;
			this->position.second = position.second;
			this->count = count;
		}
		bool operator==(const EdgePoint& rhs) 
		{ 
			return ( this->startPosition() == rhs.startPosition() and
					 this->endPosition() == rhs.endPosition() );
		}
		inline const int& startPosition() const { return position.first; }
		inline const int& endPosition() const { return position.second; }
		inline const int& edgeCount() const { return count; }
		inline const pair<int, int> positions() {return position; }
	private:
		std::pair <int, int> position;
		int count;
};
bool compareEdgePoints(const EdgePoint &lhs, const EdgePoint &rhs)
{
	return lhs.startPosition() < rhs.startPosition();
}
bool sequencesAreCloseEnough(std::string sequence1, std::string sequence2)
{//should replace this with a proper matching function but for now it checks for a max of 3 mismatches
	int mismatchLimit = 3;
	int numberOfMismatches = 0;
	#ifdef DEBUG
	std::cout << "sequencesAreCloseEnough" << std::endl << "sequence1: " << sequence1 << std::endl << "sequence2: " << sequence2 << std::endl;
	#endif
	sequence1 = boost::to_upper_copy(sequence1);
	sequence2 = boost::to_upper_copy(sequence2);
	for(int i=0; i < sequence1.size(); i++)
	{
		if(sequence1[i]!=sequence2[i])
		{
			numberOfMismatches++;
			if(numberOfMismatches > mismatchLimit)
			{	
				return false;
			}
		}
	}
	return true;
}
class SegmentsAndEdgePoints {
	public:
		SegmentsAndEdgePoints() {};
		SegmentsAndEdgePoints(stack<Segment> completedSegments, EdgePointMap completedEdgePoints)
		{
			while(!completedSegments.empty())
			{
				this->segments.push_back(Segment(completedSegments.top()));
				completedSegments.pop();
			}			
			for (EdgePointMap::iterator it = completedEdgePoints.begin(); it != completedEdgePoints.end(); ++it)
			{
				EdgePoint iteratorEdgePoint(it->first, it->second);
				this->edgePoints.push_back(EdgePoint(iteratorEdgePoint));		    	
			}
			this->sort();
		}
		SegmentsAndEdgePoints(std::string inputFileGTF, char runStrand, string chromosomeReferenceName, string test)
		{//create from GTF file which are already in SAM format		
			SegmentMap segmentMap;
			EdgePointMap edgePointMap;

			ifstream fileGTF;
			
			std::string fileSeqName, fileSource, fileFeature, fileScore, fileFrame2, fileGeneID, fileGeneIDValue, fileTranscriptID, fileTranscriptIDValue;
			char fileFrame1;

			int fileStart, fileEnd, filePreviousStart, filePreviousEnd;
			std::string fileLine;
			std::string previousFileGeneIDValue;
			int printCount = 0;
			fileGTF.open(inputFileGTF.c_str());
			printCount=0;
			
			while(getline(fileGTF, fileLine))
			{
				stringstream ss(fileLine);
				ss >> fileSeqName >> fileSource >> fileFeature >> fileStart >> fileEnd >> fileScore >> fileFrame1 >> fileFrame2 >> fileGeneID >> fileGeneIDValue >> fileTranscriptID >> fileTranscriptIDValue;
				//need to add segments and edges's to their respective vectors
				if(fileSeqName == chromosomeReferenceName and fileFeature == "exon" and fileFrame1 == runStrand)
				{
					if(fileGeneIDValue == previousFileGeneIDValue)
					{//add edge	
						EdgePoint detectedEdgePoint(make_pair(filePreviousEnd,fileStart),0);
						if(edgePointMap.find(detectedEdgePoint.positions()) == edgePointMap.end())
						{
							edgePointMap[detectedEdgePoint.positions()] = 1;
						}
						else
						{
							edgePointMap[detectedEdgePoint.positions()]++;
						}							
					}
					Segment detectedSegment(make_pair(fileStart,fileEnd));
					if(segmentMap.find(detectedSegment.positions()) == segmentMap.end())
					{
						segmentMap[detectedSegment.positions()] = 1;
					}					
					previousFileGeneIDValue = fileGeneIDValue;
					filePreviousEnd = fileEnd;				
				}
			}			
			for (SegmentMap::iterator it = segmentMap.begin(); it != segmentMap.end(); ++it)
			{
				Segment iteratorSegment(it->first);
				this->segments.push_back(Segment(iteratorSegment));		    	
			}
			for (EdgePointMap::iterator it = edgePointMap.begin(); it != edgePointMap.end(); ++it)
			{
				EdgePoint iteratorEdgePoint(it->first, it->second);
				this->edgePoints.push_back(EdgePoint(iteratorEdgePoint));		    	
			}
			fileGTF.close();		
			this->sort();
			for (std::list<Segment>::iterator segmentIterator1 = segments.begin(); segmentIterator1 != segments.end(); ++segmentIterator1)
			{
				for (std::list<Segment>::iterator segmentIterator2 = segmentIterator1; segmentIterator2 != segments.end(); ++segmentIterator2)
				{
					if(segmentIterator1==segmentIterator2) continue;
					if((*segmentIterator1).mergeIntoSegment((*segmentIterator2)))
						segmentIterator2=segments.erase(segmentIterator2);
				}
			}
		}		
		SegmentsAndEdgePoints(std::string inputFileBAM, char runStrand, string chromosomeReferenceName)
		{//create from BAM File			
			EdgePointMap edgePointMap;
			std::stack<Segment> completedSegments;
			std::stack<Segment> inProgressSegment;
			std::map<int, int> otherSegments;

			BamReader reader;
			BamAlignment alignment;	
			if (!reader.Open(inputFileBAM)) {
				std::cerr << "Could not open input BAM file." << std::endl;
				throw std::invalid_argument( "Could not open input BAM file" );
			}

			RefVector references = reader.GetReferenceData();
			int chromosomeReferenceID = chromosomeReferenceIDFromReferenceName(references, chromosomeReferenceName);
				
			while (reader.GetNextAlignment(alignment)) 
			{		
				if(runConditionsCheck(alignment,runStrand,chromosomeReferenceID)) //hard coded change later
				{		
					int startPosition = alignment.Position + 1;
					int endPosition = alignment.Position + 1;

					Segment currentCIGARSegment;
					EdgePoint detectedEdgePoint;	

					for(int i=0; i<alignment.CigarData.size(); i++) 
					{	
						int length = alignment.CigarData[i].Length;				

						switch(alignment.CigarData[i].Type)
						{					
							/*
							CIGAR allows for operations
							M 0 alignment match (can be a sequence match or mismatch)
							I 1 insertion to the reference
							D 2 deletion from the reference
							N 3 skipped region from the reference
							S 4 soft clipping (clipped sequences present in SEQ)
							H 5 hard clipping (clipped sequences NOT present in SEQ)
							P 6 padding (silent deletion from padded reference)
							= 7 sequence match
							X 8 sequence mismatch
							source: http://samtools.github.io/hts-specs/SAMv1.pdf

							PROGRAMMERS NOTE: Not all cases are considered, many handled by default case
							*/												
							case 'M'://match								
								endPosition += length - 1;
								currentCIGARSegment = Segment(make_pair(startPosition,endPosition));
								if(i==0)
								{													
									if(inProgressSegment.empty())
									{//For initialization, should only run once
										inProgressSegment.push(currentCIGARSegment);							
									}
									else if(!inProgressSegment.top().extendEndPositionWithSegment(currentCIGARSegment))
									{	
										if(!completedSegments.empty())
										{
											if(!completedSegments.top().mergeIntoSegment(inProgressSegment.top()))
											{
												completedSegments.push(inProgressSegment.top());
											}
										}
										else
										{
											completedSegments.push(inProgressSegment.top());
										}
										inProgressSegment.pop();										
										bool otherSegmentPushed = false;
										for(map<int, int>::iterator mapIterator = otherSegments.begin(); mapIterator != otherSegments.end() or mapIterator->first > currentCIGARSegment.endPosition(); ++mapIterator)
										{
											Segment iteratorSegment(std::make_pair(mapIterator->first, mapIterator->second));
											if(iteratorSegment.mergeIntoSegment(currentCIGARSegment))
											{																										
												inProgressSegment.push(iteratorSegment);
												otherSegmentPushed = true;
												otherSegments.erase(mapIterator->first);
												break;
											}
											else if(completedSegments.top().mergeIntoSegment(iteratorSegment))
											{												
												otherSegments.erase(mapIterator->first);
											}
											else if(iteratorSegment.endPosition()+1 < currentCIGARSegment.startPosition())
											{													
												completedSegments.push(Segment(make_pair(mapIterator->first, mapIterator->second)));
												otherSegments.erase(mapIterator->first);
											}
										}
										if(!otherSegmentPushed)
										{																												
											inProgressSegment.push(currentCIGARSegment);
										}
									}							
								}
								else
								{							
									bool segmentUsed = false;
									if(inProgressSegment.top().extendEndPositionWithSegment(currentCIGARSegment))
									{
										segmentUsed=true;
									}
									else
									{
										for(map<int, int>::iterator mapIterator = otherSegments.begin(); mapIterator != otherSegments.end() or mapIterator->first > currentCIGARSegment.endPosition(); ++mapIterator)
										{
											Segment iteratorSegment(std::make_pair(mapIterator->first, mapIterator->second));
											if((iteratorSegment).mergeIntoSegment(currentCIGARSegment))
											{												
												otherSegments.erase(mapIterator->first);
												otherSegments.insert(iteratorSegment.positions());
												segmentUsed=true;								
											}							
										}
									}
									if(!segmentUsed)
									{									
										otherSegments.insert(currentCIGARSegment.positions());
									}
								}						
							break;
							case 'I'://insert
								endPosition += 1;//adjustments for M length
							break;
							case 'D'://delete
								endPosition = endPosition + length + 1;//adjustments for D length and M length (+1)
							break;
							case 'N'://skipped region
								startPosition = endPosition + length + 1;
								detectedEdgePoint = EdgePoint(make_pair(endPosition, startPosition),1);
								endPosition = startPosition;
								if(edgePointMap.find(detectedEdgePoint.positions()) == edgePointMap.end())
								{
									edgePointMap[detectedEdgePoint.positions()] = 1;
								}
								else
								{
									edgePointMap[detectedEdgePoint.positions()]++;
								}						
							break;
							default:
								std::cerr << "ERROR" << std::endl;
							break;
						}
					}			
				}
			}
			reader.Close();
			//Still have a inProgressSegment and possible some in otherSegments, so add them to completed appropriately.				
			for(map<int, int>::iterator mapIterator = otherSegments.begin(); mapIterator != otherSegments.end(); ++mapIterator)
			{
				Segment iteratorSegment(std::make_pair(mapIterator->first, mapIterator->second));
				if(iteratorSegment.endPosition()+1 < inProgressSegment.top().startPosition())
				{									
					otherSegments.erase(mapIterator->first);
				}
				else if(inProgressSegment.top().mergeIntoSegment(iteratorSegment))
				{//should only occur once					
					otherSegments.erase(mapIterator->first);			
				}
				else
				{
					completedSegments.push(Segment(make_pair(mapIterator->first, mapIterator->second)));
				}								
			}
			completedSegments.push(inProgressSegment.top());
			while(!completedSegments.empty())
			{
				this->segments.push_back(Segment(completedSegments.top()));
				completedSegments.pop();
			}			
			for (EdgePointMap::iterator it = edgePointMap.begin(); it !=edgePointMap.end(); ++it)
			{
				EdgePoint iteratorEdgePoint(it->first, it->second);
				this->edgePoints.push_back(EdgePoint(iteratorEdgePoint));		    	
			}
			this->sort();
		}
		void mergeIntoSegmentsAndEdgePoints(SegmentsAndEdgePoints segmentsAndEdgePointsBeingMerged)
		{
			//So I know the list of edgepoints in each set is unique
			EdgePointMap edgePointMap;
			for (std::list<EdgePoint>::iterator edgePointIterator = this->edgePoints.begin(); edgePointIterator != this->edgePoints.end(); ++edgePointIterator)
			{				
				if(edgePointMap.find((*edgePointIterator).positions()) == edgePointMap.end())
				{
					edgePointMap[(*edgePointIterator).positions()] = (*edgePointIterator).edgeCount();
				}
				else
				{
					edgePointMap[(*edgePointIterator).positions()]+= (*edgePointIterator).edgeCount();
				}
			}
			list<EdgePoint> edgePointsBeingMerged = segmentsAndEdgePointsBeingMerged.getEdgePoints();
			for (std::list<EdgePoint>::iterator edgePointIterator = edgePointsBeingMerged.begin(); edgePointIterator != edgePointsBeingMerged.end(); ++edgePointIterator)
			{
				if(edgePointMap.find((*edgePointIterator).positions()) == edgePointMap.end())
				{
					edgePointMap[(*edgePointIterator).positions()] = (*edgePointIterator).edgeCount();
				}
				else
				{
					edgePointMap[(*edgePointIterator).positions()]+= (*edgePointIterator).edgeCount();
				}
			}
			
			list<Segment> segmentsBeingMerged = segmentsAndEdgePointsBeingMerged.getSegments();
			
			std::list<Segment>::iterator segmentIterator = this->segments.begin();
			std::list<Segment>::iterator segmentBeingMergedIterator = segmentsBeingMerged.begin();
			while(segmentIterator != this->segments.end() and segmentBeingMergedIterator != segmentsBeingMerged.end())
			{
				if((*segmentIterator).mergeIntoSegment((*segmentBeingMergedIterator)))
				{										
					segmentBeingMergedIterator=segments.erase(segmentBeingMergedIterator);					
				}
				if((*segmentIterator).endPosition() > (*segmentBeingMergedIterator).endPosition())
				{
					++segmentBeingMergedIterator;
				}
				else
				{
					++segmentIterator;
				}								
			};
			while(segmentBeingMergedIterator != segmentsBeingMerged.end())
			{
				this->segments.push_back((*segmentBeingMergedIterator));
				++segmentBeingMergedIterator;
			}

			this->edgePoints.clear();
			for (EdgePointMap::iterator it = edgePointMap.begin(); it != edgePointMap.end(); ++it)
			{
				EdgePoint iteratorEdgePoint(it->first, it->second);
				this->edgePoints.push_back(EdgePoint(iteratorEdgePoint));		    	
			}
		}
		void sort() 
		{
			this->segments.sort(compareSegments);
			this->edgePoints.sort(compareEdgePoints);
		}
		std::list<Segment> getSegments() const { return segments; }
		std::list<EdgePoint> getEdgePoints() const { return edgePoints; }
		void printSegments()
		{							
			std::cout << std::endl << "*** Segments ***" << std::endl;
			int numberOfSegments = 0;
			for (std::list<Segment>::iterator segmentIterator = segments.begin(); segmentIterator != segments.end(); ++segmentIterator)
			{
				numberOfSegments++;
				std::cout << numberOfSegments << ") " << (*segmentIterator).startPosition() << "," << (*segmentIterator).endPosition() << '\t';
				if(numberOfSegments%4==0)
				{
					std::cout << std::endl;
				}
			}		
			std::cout << std::endl << "Total number of Segments: " << numberOfSegments << std::endl;	
		}
		void printEdgePoints()
		{			
			std::cout <<  std::endl << "*** EdgePoints ***" << std::endl;
			int numberOfEdgePoints = 0;
				for (std::list<EdgePoint>::iterator edgePointIterator = edgePoints.begin(); edgePointIterator != edgePoints.end(); ++edgePointIterator)
				{
					numberOfEdgePoints++;
					std::cout << numberOfEdgePoints << ") "<< (*edgePointIterator).startPosition() << "->" << (*edgePointIterator).endPosition() << " EC:" << (*edgePointIterator).edgeCount() << '\t';
					if(numberOfEdgePoints%3==0)
					{
						std::cout << std::endl;
					}
				}
			std::cout << std::endl << "Total number of EdgePoints: " << numberOfEdgePoints << std::endl;	
		}
	private:
		std::list<Segment> segments;
		std::list<EdgePoint> edgePoints;		
};
class ReferenceFASTAData
{
	public:
		ReferenceFASTAData(string referenceFileFASTA)
		{
			referenceFASTAData = getFileContents(referenceFileFASTA.c_str());//contains reference
			referenceFASTAData.erase(referenceFASTAData.begin(),referenceFASTAData.begin()+referenceFASTAData.find("\n"));//remove header line
			referenceFASTAData.erase(std::remove_if(referenceFASTAData.begin(),referenceFASTAData.end(), ::isspace ), referenceFASTAData.end());//remove whitespace
		}
		std::string obtainSubSequence(int startPosition, int endPosition)
		{//stored in base 0 while positions are in base 1
			return referenceFASTAData.substr(startPosition - 1, endPosition - startPosition + 1);
		}		
	private:
		std::string referenceFASTAData;
};
class SpliceGraph
{//apply smart logic here	
	public:
		SpliceGraph(SegmentsAndEdgePoints segmentsAndEdgePoints)
		{			
			std::list<Segment> segments = segmentsAndEdgePoints.getSegments();
			std::list<EdgePoint> edgePoints = segmentsAndEdgePoints.getEdgePoints();
			
			for (std::list<Segment>::iterator segmentIterator = segments.begin(); segmentIterator != segments.end(); ++segmentIterator)
			{
				for (std::list<EdgePoint>::iterator edgePointIterator = edgePoints.begin(); edgePointIterator != edgePoints.end(); ++edgePointIterator)
				{
					if(	(*segmentIterator).startPosition() < (*edgePointIterator).startPosition() and
						(*segmentIterator).endPosition() > (*edgePointIterator).startPosition() )
					{//123|456 means | is at 3 so 123 and 456 												
						segments.push_back(Segment(make_pair((*edgePointIterator).startPosition()+1,(*segmentIterator).endPosition())));
						(*segmentIterator) = Segment(make_pair((*segmentIterator).startPosition(),(*edgePointIterator).startPosition()));
						EdgePoint newEdgePoint(make_pair((*edgePointIterator).startPosition(),(*edgePointIterator).startPosition()+1),0);						
						edgePoints.push_back(newEdgePoint);
					}
					if(	(*segmentIterator).startPosition() < (*edgePointIterator).endPosition() and
						(*segmentIterator).endPosition() > (*edgePointIterator).endPosition())
					{//789|0123 means | is at 0 so 789 and 0123
						segments.push_back(Segment(make_pair((*edgePointIterator).endPosition(),(*segmentIterator).endPosition())));
						(*segmentIterator) = Segment(make_pair((*segmentIterator).startPosition(),(*edgePointIterator).endPosition()-1));
						EdgePoint newEdgePoint(make_pair((*edgePointIterator).endPosition()-1,(*edgePointIterator).endPosition()),0);
						edgePoints.push_back(newEdgePoint);
					}
				}
			}
			for (std::list<Segment>::iterator segmentIterator = segments.begin(); segmentIterator != segments.end(); ++segmentIterator)
			{
				spliceGraphVertex newVertex = add_vertex(spliceGraph);
				std::stringstream name;
				name << "chr3 | |" << (*segmentIterator).startPosition()  << "-" << (*segmentIterator).endPosition();//+1 is to change from = start BAm format to 1 start SAM format				
				spliceGraph[newVertex].startPosition = (*segmentIterator).startPosition();
				spliceGraph[newVertex].endPosition = (*segmentIterator).endPosition();
				spliceGraph[newVertex].name = name.str();				
			}
			#ifdef DEBUG
			printCurrentRunTime(startTime, processTime, totalTime);	
			#endif
			for (std::list<EdgePoint>::iterator edgePointIterator = edgePoints.begin(); edgePointIterator != edgePoints.end(); ++edgePointIterator)
			{		
				graph_traits <SpliceGraphType>::vertex_iterator vi, vi_end;
				graph_traits <SpliceGraphType>::vertex_iterator vi2, vi2_end;
				for (boost::tie(vi, vi_end) = vertices(spliceGraph); vi != vi_end; ++vi) 					
				{					
					if(spliceGraph[*vi].endPosition == (*edgePointIterator).startPosition())
					{
						break;
					}
				}
				for (boost::tie(vi2, vi2_end) = vertices(spliceGraph); vi2 != vi2_end; ++vi2) 
				{
					if(spliceGraph[*vi2].startPosition == (*edgePointIterator).endPosition())
					{
						break;
					}
				}				
				add_edge(*vi,*vi2,spliceGraph);					
			}
			#ifdef DEBUG
			printCurrentRunTime(startTime, processTime, totalTime);	
			#endif
		}		
		int spliceGraphRemoveMappedVerticesErrors(std::string inputFileBAM, char runStrand, string chromosomeReferenceName, ReferenceFASTAData referenceFASTAData) 
		{
			std::list<spliceGraphVertex> segmentsToCheckPlus;
			std::list<spliceGraphVertex> segmentsToCheckMinus;

			BamReader reader;
			BamAlignment alignment;

			if (!reader.Open(inputFileBAM)) {
				std::cerr << "Could not open input BAM file." << std::endl;
				return 1;
			}

			RefVector references = reader.GetReferenceData();
			int chromosomeReferenceID = chromosomeReferenceIDFromReferenceName(references, chromosomeReferenceName);

			/*
			Mark all edges that need to be checked. The criteria for being checked is that current node has to be 1 '
			character away from another node and that it only has 1 incoming or 1 outgoing edge.
			*/	
			graph_traits <SpliceGraphType>::vertex_iterator vi, vi_end;
			graph_traits <SpliceGraphType>::edge_iterator ei, ei_end;	
			graph_traits <SpliceGraphType>::out_edge_iterator oei, oei_end;
			graph_traits <SpliceGraphType>::in_edge_iterator iei, iei_end;

			for(boost::tie(ei, ei_end) = edges(spliceGraph); ei != ei_end; ++ei) 
			{
				if(spliceGraph[source(*ei, spliceGraph)].endPosition+1==spliceGraph[target(*ei, spliceGraph)].startPosition) 
				{
					//Check for target or +edge	
					int outgoingEdges = 0;
					int incomingEdges = 0;			
					for (boost::tie(oei, oei_end) = out_edges(target(*ei, spliceGraph), spliceGraph); oei != oei_end; ++oei)
					{
						outgoingEdges++;
						if(outgoingEdges>0)
						{
							break;
						}
					}
					for (boost::tie(iei, iei_end) = in_edges(target(*ei, spliceGraph), spliceGraph); iei != iei_end; ++iei)
					{
						incomingEdges++;
						if(incomingEdges>1)
						{
							break;
						}
					}
					if(outgoingEdges==0 and incomingEdges==1)
					{						
						segmentsToCheckPlus.push_back(target(*ei, spliceGraph));
					}
				}
				if(spliceGraph[source(*ei, spliceGraph)].endPosition+1==spliceGraph[target(*ei, spliceGraph)].startPosition) 
				{		
					int outgoingEdges = 0;
					int incomingEdges = 0;			
					for (boost::tie(oei, oei_end) = out_edges(source(*ei, spliceGraph), spliceGraph); oei != oei_end; ++oei)
					{
						outgoingEdges++;
						if(outgoingEdges>1)
						{
							break;
						}
					}
					for (boost::tie(iei, iei_end) = in_edges(source(*ei, spliceGraph), spliceGraph); iei != iei_end; ++iei)
					{
						incomingEdges++;
						if(incomingEdges>0)
						{
							break;
						}
					}
					if(outgoingEdges==1 and incomingEdges==0)
					{					
						segmentsToCheckMinus.push_back(source(*ei, spliceGraph));
					}
				}
			}			
			segmentsToCheckPlus.sort(boost::bind(&SpliceGraph::compareStartPositions, this, _1, _2));//bind in order class to indentify function
			segmentsToCheckMinus.sort(boost::bind(&SpliceGraph::compareStartPositions, this, _1, _2));//bind in order class to indentify function
			
			cout << "PLUS" << endl;
			for (std::list<spliceGraphVertex>::iterator segmentIterator = segmentsToCheckPlus.begin(); segmentIterator != segmentsToCheckPlus.end(); ++segmentIterator)
			{
				cout << spliceGraph[(*segmentIterator)].name << endl;				
			}
			cout << "MINUS" << endl;
			for (std::list<spliceGraphVertex>::iterator segmentIterator = segmentsToCheckMinus.begin(); segmentIterator != segmentsToCheckMinus.end(); ++segmentIterator)
			{
				cout << spliceGraph[(*segmentIterator)].name << endl;				
			}
			
			std::list<spliceGraphVertex>::iterator segmentIteratorPlus = segmentsToCheckPlus.begin();
			std::list<spliceGraphVertex>::iterator segmentIteratorMinus = segmentsToCheckMinus.begin();			
			
			while (reader.GetNextAlignment(alignment) and segmentIteratorPlus != segmentsToCheckPlus.end() and segmentIteratorMinus != segmentsToCheckMinus.end()) 
			{
				if(runConditionsCheck(alignment,runStrand,chromosomeReferenceID))
				{
					int startPosition = alignment.Position + 1;
					int endPosition = alignment.Position + 1;					

					if(alignment.Position+1 > spliceGraph[(*segmentIteratorPlus)].endPosition)
					{											
						//segmentsToCheckPlus.erase(segmentIteratorPlus);
						//segmentIteratorPlus = segmentsToCheckPlus.begin();
					}
					if(alignment.Position+1 > spliceGraph[(*segmentIteratorMinus)].endPosition)
					{												
						//segmentIteratorMinus = segmentsToCheckMinus.erase(segmentIteratorMinus);
						//segmentIteratorMinus = segmentsToCheckMinus.begin();
					}
									
					for(int i=0; i<alignment.CigarData.size(); i++) 
					{				
						int length = alignment.CigarData[i].Length;
						switch(alignment.CigarData[i].Type)
						{
							case 'M'://match								
								endPosition += length - 1;
								if(segmentIteratorPlus != segmentsToCheckPlus.end() and endPosition > spliceGraph[(*segmentIteratorPlus)].startPosition)
								{//Always at the end of a fragment																		
									for(std::list<spliceGraphVertex>::iterator segmentIterator = segmentsToCheckPlus.begin(); segmentIterator != segmentsToCheckPlus.end() and spliceGraph[(*segmentIterator)].startPosition < endPosition; ++segmentIterator)
									{
										int segmentLength = spliceGraph[(*segmentIterator)].endPosition - spliceGraph[(*segmentIterator)].startPosition + 1;
										if(	spliceGraph[(*segmentIterator)].endPosition == endPosition and 
											alignment.Length > segmentLength )
										{//This check should obtain all plus checks											
											int segmentPosition = 0;//to track position on readSequence
											std::string readSequence = alignment.QueryBases;											
											//Don't want to store all reads, just those affected so need to loop through CIGAR again)
											for(int j=0; j<alignment.CigarData.size(); j++)
											{
												switch(alignment.CigarData[j].Type)
												{
													case 'M':
														segmentPosition += alignment.CigarData[j].Length;
														if(i==j)
														{
															if((segmentPosition - segmentLength + 1) < readSequence.size())
															{//error check to ensure substr works
																std::string compareSequence = readSequence.substr(segmentPosition - segmentLength, segmentLength);
																std::string referenceSequence = referenceFASTAData.obtainSubSequence(spliceGraph[(*segmentIterator)].startPosition,spliceGraph[(*segmentIterator)].endPosition);
																if(	spliceGraph[(*segmentIterator)].name == "REMOVE")
																{
																	break;
																}
																else if(compareSequence.length() > 0 and 
																		compareSequence.length() == referenceSequence.length() and 
																		sequencesAreCloseEnough(compareSequence,referenceSequence) )
																{																
																	spliceGraph[(*segmentIterator)].name = "REMOVE";
																	/*segmentsToCheckPlus.erase(segmentIterator);
																	segmentIterator = segmentsToCheckPlus.begin();
																	if(spliceGraph[(*segmentIteratorPlus)].name == "REMOVE")
																	{
																		segmentIteratorPlus = segmentsToCheckPlus.begin();	
																	}*/
																	break;
																}
															}
														}
													break;
													case 'I':														
														readSequence.erase(segmentPosition,alignment.CigarData[j].Length);														
													break;
													case 'D':
														readSequence.insert(segmentPosition,alignment.CigarData[j].Length,'*');
														segmentPosition += alignment.CigarData[j].Length;
													break;
													case 'N':
														readSequence.erase(0,segmentPosition);
														segmentPosition = 0;
													break;
												}
											}
										}
									}
								}

								if(segmentIteratorMinus != segmentsToCheckMinus.end() and startPosition >= spliceGraph[(*segmentIteratorMinus)].startPosition)
								{//always at the start of a fragment
									for(std::list<spliceGraphVertex>::iterator segmentIterator = segmentsToCheckMinus.begin(); segmentIterator != segmentsToCheckMinus.end() and spliceGraph[(*segmentIterator)].startPosition < endPosition; ++segmentIterator)
									{
										int segmentLength = spliceGraph[(*segmentIterator)].endPosition - spliceGraph[(*segmentIterator)].startPosition + 1;
										if(	spliceGraph[(*segmentIterator)].startPosition == startPosition and 
											alignment.Length > segmentLength )
										{//This check should obtain all minus checks											
											int segmentPosition = 0;//to track position on readSequence
											std::string readSequence = alignment.QueryBases;											

											//Don't want to store all reads, just those affected so need to loop through CIGAR again)
											for(int j=0; j<alignment.CigarData.size(); j++)
											{
												switch(alignment.CigarData[j].Type)
												{
													case 'M':																										
														segmentPosition +=alignment.CigarData[j].Length;
														if(spliceGraph[(*segmentIterator)].startPosition+segmentPosition > spliceGraph[(*segmentIterator)].endPosition)
														{
																
															if(segmentLength < readSequence.size())
															{//error check to ensure substr works
																std::string compareSequence = readSequence.substr(0, segmentLength);
																std::string referenceSequence = referenceFASTAData.obtainSubSequence(spliceGraph[(*segmentIterator)].startPosition,spliceGraph[(*segmentIterator)].endPosition);
																if(	spliceGraph[(*segmentIterator)].name == "REMOVE")
																{
																	break;
																}
																else if(compareSequence.length() > 0 and 
																		compareSequence.length() == referenceSequence.length() and 
																		sequencesAreCloseEnough(compareSequence,referenceSequence) )
																{															
																	spliceGraph[(*segmentIterator)].name = "REMOVE";
																	
																/*	segmentsToCheckMinus.erase(segmentIterator);
																	segmentIterator = segmentsToCheckMinus.begin();
																	if(spliceGraph[(*segmentIteratorMinus)].name=="REMOVE")
																	{
																		segmentIteratorMinus = segmentsToCheckMinus.begin();
																	}*/
																	break;
																}
															}
														}														
													break;
													case 'I':
														readSequence.erase(segmentPosition,alignment.CigarData[j].Length);														
													break;
													case 'D':
														readSequence.insert(segmentPosition,alignment.CigarData[j].Length,'*');
														segmentPosition += alignment.CigarData[j].Length;
													break;
													case 'N':
														readSequence.erase(0,segmentPosition);
														segmentPosition = 0;
													break;
												}
											}
										}
									}
								}
							break;
							case 'I'://insert
								endPosition += 1;//adjustments for M length
							break;
							case 'D'://delete
								endPosition = endPosition + length + 1;//adjustments for D length and M length (+1)
							break;
							case 'N'://skipped region
								startPosition = endPosition + length + 1;
								endPosition = startPosition;
							break;
							default:
								std::cerr << "ERROR" << std::endl;
							break;
						}	
					}
				}
			}
			reader.Close();			

			return 0;	
		} 
		void benchmarkAgainstGTF(std::string benchmarkFile, char runStrand, string chromosomeReferenceName)
		{//The GTF file is assumed to be ordered such genes are grouped together in an order of increasing position (earlier positions being in higher rows)					

			vector <int> gtfEdges;//shouldn't be needed as a global now as I should be able to merge graphs for the check
			ifstream fileGTF;
			std::string fileSeqName, fileSource, fileFeature, fileScore, fileFrame2, fileGeneID, fileGeneIDValue, fileTranscriptID, fileTranscriptIDValue;
			int fileStart, fileEnd, filePreviousStart, filePreviousEnd;
			char fileFrame1;
			std::string fileLine;
			std::string previousFileGeneIDValue;
			int printCount = 0;
			fileGTF.open(benchmarkFile.c_str());
			printCount=0;
			
			while(getline(fileGTF, fileLine))
			{
				stringstream ss(fileLine);
				ss >> fileSeqName >> fileSource >> fileFeature >> fileStart >> fileEnd >> fileScore >> fileFrame1 >> fileFrame2 >> fileGeneID >> fileGeneIDValue >> fileTranscriptID >> fileTranscriptIDValue;
				if(fileSeqName == chromosomeReferenceName and fileFeature == "exon" and runStrand == fileFrame1)
				{
					if(fileGeneIDValue == previousFileGeneIDValue)
					{	
						//Print connections			
						std::cout << printCount << ")\t" << filePreviousEnd << "->" << fileStart << "\t";
						if(printCount%4==3)
							std::cout << std::endl;

						bool edgeFound = false;
						graph_traits <SpliceGraphType>::edge_iterator ei, ei_end;
						for (boost::tie(ei, ei_end) = edges(spliceGraph); ei != ei_end; ++ei) 
						{
							if(	spliceGraph[source(*ei, spliceGraph)].endPosition == filePreviousEnd or
								spliceGraph[target(*ei, spliceGraph)].startPosition == fileStart )
								{
									edgeFound=true;
								}
						}

						if(!edgeFound)
						{				
							gtfEdges.push_back(0);
						}
						else
						{			
							gtfEdges.push_back(1);
						}
						filePreviousEnd = fileEnd;			
						printCount++;			
					}
					previousFileGeneIDValue = fileGeneIDValue;
					filePreviousEnd = fileEnd;
					filePreviousStart = fileStart;
				}
			}
			fileGTF.close();

			std::cout << std::endl;
			for(int i=0; i<gtfEdges.size(); i++)
			{//Display results
				if(gtfEdges[i] == 1)
					std::cout << i << ")\t" << gtfEdges[i] << "\t";
				else
					std::cout << i << ")\t" << gtfEdges[i] << " <--\t";
				if(i%7==6)
					std::cout << std::endl;
			}
			std::cout << std::endl;

		}
		void printGraph()
		{
			boost::print_graph(spliceGraph, get(&MyVertex::name,spliceGraph));
		}	
		void generateDOTFile(string outputFileName)
		{
			ofstream outstream;		
			outstream.open(outputFileName.c_str());
			boost::write_graphviz(outstream, spliceGraph, boost::make_label_writer(boost::get(&MyVertex::name,spliceGraph)));
			outstream.close();
			std::cout << outputFileName << " made" << std::endl;
		}		
	private:
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, MyVertex, MyEdge> SpliceGraphType;
		typedef boost::graph_traits<SpliceGraphType>::vertex_descriptor spliceGraphVertex;
		typedef boost::graph_traits<SpliceGraphType>::edge_descriptor spliceGraphEdge;
		SpliceGraphType spliceGraph;					
		bool compareStartPositions(spliceGraphVertex lhs, spliceGraphVertex rhs)
		{
			return spliceGraph[lhs].startPosition < spliceGraph[rhs].endPosition;
		}
};

//Debug Function
int statsOnInputBAMFile(std::string inputFileBAM)
{
	#ifdef DEBUG
	BamReader reader;
	BamAlignment alignment;	
	if (!reader.Open(inputFileBAM)) {
		std::cerr << "Could not open input BAM file." << std::endl;
		return 1;
	}

	int numberOfAlignments=0;
	int numberOfUsedAlignments=0;
	int numberOfPlusStrandCIGAR_M=0;
	int numberOfPlusStrandCIGAR_N=0;
	int numberOfMinusStrandCIGAR_M=0;
	int numberOfMinusStrandCIGAR_N=0;

	while (reader.GetNextAlignment(alignment)) 
	{
		numberOfAlignments++;
		if(runConditionsCheck(alignment,'o',-1))
		{
			numberOfUsedAlignments++;
			for(int i=0; i<alignment.CigarData.size(); i++) 
			{
				switch(alignment.CigarData[i].Type)
				{										
					case 'M'://match
						if(runConditionsCheck(alignment,'+',-1))
							numberOfPlusStrandCIGAR_M++;		
						if(runConditionsCheck(alignment,'-',-1))
							numberOfMinusStrandCIGAR_M++;
					break;
					case 'N':
						if(runConditionsCheck(alignment,'+',-1))
							numberOfPlusStrandCIGAR_N++;
						if(runConditionsCheck(alignment,'-',-1))
							numberOfMinusStrandCIGAR_N++;
					break;
				}
			}
		}
	}

	std::cout 	<< " Number of Alignments :\t" << numberOfAlignments << std::endl
				<< " Number of CIGAR reads:\t" << numberOfPlusStrandCIGAR_M + numberOfMinusStrandCIGAR_M << std::endl
				<< "          +CIGAR reads:\t" << numberOfPlusStrandCIGAR_M << std::endl
				<< "          -CIGAR reads:\t" << numberOfMinusStrandCIGAR_M << std::endl				
				<< " Number of Edge Points:\t" << numberOfPlusStrandCIGAR_N + numberOfMinusStrandCIGAR_N << std::endl
				<< "          +Edge Points:\t" << numberOfPlusStrandCIGAR_N << std::endl
				<< "          -Edge Points:\t" << numberOfMinusStrandCIGAR_N << std::endl
				<< std::endl;
	reader.Close();
	#endif
	return 0;
}

//Debug Function: to check if a read covers a position
int findReadsCoveringPosition(std::string inputFileBAM, int positionToCheckForReads) 
{
	#ifdef DEBUG
	BamReader reader;
	BamAlignment alignment;
	
	if (!reader.Open(inputFileBAM)) {
		std::cerr << "Could not open input BAM file." << std::endl;
		return 1;
	}

	positionToCheckForReads -= 1;//as this number to check is Base 1 and the BAM file is base 0
	std::cout << "*** Reads covering position: " << positionToCheckForReads+1 << " ***" << std::endl;
	while (reader.GetNextAlignment(alignment)) 
	{
		int startPosition = alignment.Position;
		int endPosition = alignment.Position;
		std::string strand;
		if(	(alignment.IsFirstMate() and alignment.IsReverseStrand()) or 
			(!alignment.IsFirstMate() and !alignment.IsReverseStrand())) 
		{
			strand = '+';
		} 
		else if((alignment.IsFirstMate() and !alignment.IsReverseStrand()) or 
				(!alignment.IsFirstMate() and alignment.IsReverseStrand()))  
		{
			strand = '-';
		}
		if(runConditionsCheck(alignment,'o',-1))//conditions check forward and back
		{
			for(int i=0; i<alignment.CigarData.size(); i++) 
			{
				switch(alignment.CigarData[i].Type)
				{										
					case 'M'://match
						endPosition += alignment.CigarData[i].Length - 1;
						if(startPosition <= positionToCheckForReads and endPosition >= positionToCheckForReads)
						{								
								std::cout 	<< alignment.Name <<"\t: " << alignment.QueryBases << std::endl;
								std::cout 	<< strand << "CIGAR\t\t: ";
								for(int j=0; j<alignment.CigarData.size(); j++) 
								{
									std::cout << alignment.CigarData[j].Length << alignment.CigarData[j].Type;
								}
								std::cout << "\tAlignment Start: " << alignment.Position+1 << "\tSegment Start: " << startPosition+1 
									 << " End: " << endPosition+1 << "\t" << alignment.CigarData[i].Length;
								std::cout << std::endl;
						}
					break;
					case 'I'://insert
					break;
					case 'D'://delete
						endPosition += alignment.CigarData[i].Length;//Adjust position as provided by cigar (relative position)				
					break;
					case 'N'://skipped region
						startPosition = endPosition + alignment.CigarData[i].Length + 1;								
						endPosition = startPosition;
					break;
					default:
						std::cerr << "ERROR" << std::endl;
					break;
				}	
			}
		}		
	}
	std::cout << std::endl;
	reader.Close();
	#endif
	return 0;
}

int main(int argc, char* argv[])
{
	std::cout << "****** Graphinator start ******" << std::endl;

	#ifdef DEBUG	
	time(&startTime);
	time(&processTime);
	#endif

	std::string inputFile;
	if ( argc == 2 ) 
	{
		inputFile = argv[1] ;
	}
	else 
	{
		std::cerr << "Wrong number of arguments." << std::endl;
		return 1;
	}
	
	string inputFileBAM = "/home/kimng/Code/Graphanator/build/accepted_hits_larger.bam";//"/demeter-internal/jasi/projects/splice_graph/chr3_perfect.bam";
	string inputFileFASTA = "/demeter-internal/jasi/projects/splice_graph/chr3.fa";
	string inputFileGTF = "/home/kimng/Code/Graphanator/build/RefSeq_chr3_23M_26M.gtf";//"sftp://kimng@demeter/demeter-internal/jasi/projects/splice_graph/chr3.gtf";
	char runStrand = '+';
	string chromosomeReferenceName = "chr3";
	string outputFileSpliceGraphDOT = "test_graph_plus.dot";
	string outputFileSpliceGraphPDF = "test_graph_plus.pdf";
	string outputFileFilteredGraphDOT = "test_graph_plus_filtered.dot";
	string outputFileFilteredGraphPDF = "test_graph_plus_filtered.pdf";

	statsOnInputBAMFile(inputFileBAM);
	printCurrentRunTime(startTime, processTime, totalTime);	
	SegmentsAndEdgePoints segmentsAndEdgePoints;
	SegmentsAndEdgePoints segmentsAndEdgePointsFromGTFPlus(inputFileGTF, '+', chromosomeReferenceName, "junk"); 	
	//SegmentsAndEdgePoints segmentsAndEdgePointsFromGTFMinus(inputFileGTF, '-', chromosomeReferenceName, "junk"); 
	//SegmentsAndEdgePoints segmentsAndEdgePointsFromBAM(inputFileBAM, runStrand, chromosomeReferenceName);
	printCurrentRunTime(startTime, processTime, totalTime, "- SegmentsAndEdgePoints made from BAM");
	segmentsAndEdgePoints.mergeIntoSegmentsAndEdgePoints(segmentsAndEdgePointsFromGTFPlus);
	
	//segmentsAndEdgePoints.mergeIntoSegmentsAndEdgePoints(segmentsAndEdgePointsFromGTFMinus);//Test out merging the 2 GTF's (each strand), was doing testing on minus where it works and having issues with plus
	//segmentsAndEdgePoints.mergeIntoSegmentsAndEdgePoints(segmentsAndEdgePointsFromBAM);//somethings off with pointer logic so documented out for now
	segmentsAndEdgePoints.printSegments();
	segmentsAndEdgePoints.printEdgePoints();
	printCurrentRunTime(startTime, processTime, totalTime);
	SpliceGraph spliceGraph(segmentsAndEdgePoints);	
	printCurrentRunTime(startTime, processTime, totalTime, "- SegmentsAndEdgePoints converted to SpliceGraph");
	ReferenceFASTAData referenceFASTAData(inputFileFASTA);
	printCurrentRunTime(startTime, processTime, totalTime, "- ReferenceFASTAData put into Memory");
	spliceGraph.spliceGraphRemoveMappedVerticesErrors(inputFileBAM, runStrand, chromosomeReferenceName, referenceFASTAData);
	printCurrentRunTime(startTime, processTime, totalTime, "- spliceGraphRemoveMappedVerticesErrors complete");
	spliceGraph.generateDOTFile("test.dot");
	//spliceGraph.printGraph();
	printCurrentRunTime(startTime, processTime, totalTime, "- dot file made");
	//spliceGraph.benchmarkAgainstGTF(inputFileGTF, runStrand, chromosomeReferenceName);
	//printCurrentRunTime(startTime, processTime, totalTime);
	string systemCommandToRun = "dot -Tpdf -o test.pdf test.dot";
	system(systemCommandToRun.c_str());//run on system
	printCurrentRunTime(startTime, processTime, totalTime);

	//findReadsCoveringPosition(inputFileBAM, 25835578);

	std::cout << "****** Graphinator end   ******" << std::endl;
	return 0;
}
