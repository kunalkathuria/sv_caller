#Identify all discordant reads
#$ Change "tid" to "ref_name"
import math
import sys
import pysam
import time
from itertools import izip

PERMUTATION_THRESH = int(sys.argv[3])#20
CALC_THRESH = int(sys.argv[6]) #500000
MATCH_PCT = int(sys.argv[7]) #0
MATCH_THRESH = float(sys.argv[5]) #.9 # RELATIVE ALIGNMENT PRECISION OF PRIMARY VS SECONDARY ALIGNMENT
PCT_THRESH = float(sys.argv[4]) #1
FILE1 = "../data/bams/aln1s.bam" # read 1 discordants
FILE2 = "../data/bams/aln2s.bam"# read 2 discordants
BAMFILE1 = sys.argv[1] #name-sorted BAM file. Make sure index file has same name with ".bai" suffix added.
BAMFILE2 = sys.argv[2] # position-sorted same BAM file.

def calcMeanSig(file1, file2):
   bamfile = pysam.Samfile(file1,"rb")
   summedIL = 0
   summedQL = 0
   counter = 0
   stdev = 0
   max_IL = 0

   while counter < CALC_THRESH:

        try:

            m1 = bamfile.next()
            m2 = bamfile.next()

        except StopIteration:
            break

        if m1.query_name == m2.query_name and m1.is_proper_pair and m2.is_proper_pair and not m1.is_secondary and not m1.is_supplementary and not m2.is_secondary and not m2.is_supplementary:

             if m1.reference_start < m2.reference_start:
                 minm = m1
                 maxm = m2
             else:
                 maxm = m1
                 minm = m2

	     IL = maxm.reference_end - minm.reference_start	
	     if IL > max_IL:
		max_IL = IL
             summedIL+= IL
	     summedQL+= m1.infer_query_length() + m2.infer_query_length()
             counter+= 1

   meanIL = summedIL/counter
   meanQL = summedQL/(2*counter)
   bamfile.close()
   bamfile = pysam.Samfile(file1,"rb")

   counter=0
   while counter < CALC_THRESH:

        try:

            m1 = bamfile.next()
            m2 = bamfile.next()

        except StopIteration:
            break

        if m1.query_name == m2.query_name and m1.is_proper_pair and m2.is_proper_pair and not m1.is_secondary and not m1.is_supplementary and not m2.is_secondary and not m2.is_supplementary:

             if m1.reference_start < m2.reference_start:
                 minm = m1
                 maxm = m2
             else:
                 maxm = m1
                 minm = m2

             stdev+= (maxm.reference_end - minm.reference_start - meanIL)**2
	     counter+=1

   bamfile.close()
   bamfile = pysam.AlignmentFile(file2,"rb")

   cov = 0
   width = 20
   counter2=0

   for pileupcolumn in bamfile.pileup():

   	cov+=pileupcolumn.n
	counter2+=1
	if counter2 > CALC_THRESH:
		break

   bamfile.close()
   stdev = stdev/counter
   stdev = stdev**(.5)
   cov=cov/counter2

   print meanQL, meanIL, stdev, cov, max_IL
   fp = open("../results/text/bam_stats.txt","w")
   fp.write("%s\n%s\n%s\n%s\n%s\n" %(meanQL, meanIL, stdev,cov, max_IL))
   fp.close()
   return meanQL, meanIL, stdev

[RDL, MEAN_D, SIG_D] = calcMeanSig(BAMFILE1, BAMFILE2)

class Cluster(object):

    def __init__(self):

       #Whether l_bound is start or end of left-aligned read (in reference) depends on orientation of read in reference
       self.l_bound=None
       self.r_bound=None
       self.l_bound_O=None # debug variable. first pair in cluster: all positions assessed relative to this one for future additions to cluster
       self.r_bound_O=None
       self.C_type=None # LR orientation: F=0, R=1
       self.Cmap_type = None # 0 if both reads map, 1 if one read maps only
       self.ltid = None
       self.rtid = None
       self.l_target = None
       self.r_target = None

    def __str__(self):
        return "%s %s %s %s %s %s %s" % (self.ltid, self.l_bound, self.rtid, self.r_bound, self.C_type, self.l_target, self.r_target)


def findNumberMatches(cigartups):
    nummatches = 0
    if cigartups == None:
        return 0
    for op,oplen in cigartups:
        if op in [0,1,7]:
            nummatches += int(oplen)

    return nummatches

def FormDiscordant(list1, list2, DList1, DList2):

    counter = 0
    OUTER_D = MEAN_D # non-discordant default value. safe.

    first = 1
    
    for al1 in list1:

            for al2 in list2:
		
                counter = counter + 1
                al1_match = 0
		al2_match = 0
		map_type = 3 # improper value
                al1_reverse = 2 # improper value
                al2_reverse = 2 # improper value
                CLType = "22" # "(map_type)(L is reverse?)(R is reverse?)(R > L pos?)" # improper value
                indivPos = -1 # temp position variable for map-type 1 clusters
                lf_bound = -1
                rt_bound = -1
                l_orient = 2
                r_orient = 2
                indivOrient = 2 
		lf_target = -1
		rt_target = -1
	
		if counter > PERMUTATION_THRESH:
			DList1 = []
			DList2 = []
			return

                if first == 1:
                    MatchF1 = findNumberMatches(al1.cigartuples)
                    MatchF2 = findNumberMatches(al2.cigartuples)
                    ASF1 = int(al1.get_tag("AS"))
                    ASF2 = int(al2.get_tag("AS"))
                    

                NMatch1 = findNumberMatches(al1.cigartuples)
                NMatch2 = findNumberMatches(al2.cigartuples)

                if NMatch1 > MATCH_PCT*RDL and float(NMatch1)/float(MatchF1) >= MATCH_THRESH and float(al1.get_tag("AS")) >= PCT_THRESH*ASF1:
                    al1_match = 1

                else:
                    continue # can use return here if sure BAM file reports alignments in order of decreasing AS
		
                if NMatch2 > MATCH_PCT*RDL and float(NMatch2)/float(MatchF2) >= MATCH_THRESH and float(al2.get_tag("AS")) >= PCT_THRESH*ASF2:
                    al2_match = 1

                else:
                    continue
		
                if first == 1:
                    first = 0

                # al_match check in conditionals is now redundant but safe in case needed later
                
                if (al1.reference_start== None or al1.reference_start <= -1 or al1.is_unmapped) and (al2.reference_start!=None and al2.reference_start > -1 and al2_match):

                    map_type = 1
                    
                    # left tid set by convention for one mate mappings
	            left_tid = al2.reference_name

                    if al2.is_reverse:
 
                        indivPos = al2.reference_start
                        indivOrient = 1
                        OUTER_D = -1

			# based on likely situation for random insertion
                        lf_target = al2.query_alignment_start
                        rt_target = al1.query_alignment_end
                        
                    else:

                        indivPos = al2.reference_end
                        indivOrient = 0
                        OUTER_D = -1

			# based on likely situation for random insertion
                        lf_target = al2.query_alignment_end
                        rt_target = al1.query_alignment_start
                            

                if (al2.reference_start== None or al2.reference_start <= -1 or al2.is_unmapped) and (al1.reference_start!=None and al1.reference_start > -1 and al1_match):

                    map_type = 1

		    left_tid = al1.reference_name	
	
                    if al1.is_reverse:
 
                        indivPos = al1.reference_start
                        indivOrient = 1
                        OUTER_D = -1

			# based on likely situation for random insertion
			lf_target = al1.query_alignment_start
			rt_target = al2.query_alignment_end
                        
                    else:

                        indivPos = al1.reference_end
                        indivOrient = 0
                        OUTER_D = -1

			# based on likely situation for random insertion
                        lf_target = al1.query_alignment_end
                        rt_target = al2.query_alignment_start
  
                if al1.reference_start!= None and al1.reference_start > -1 and al2.reference_start!=None and al2.reference_start > -1 and al1_match and al2_match:
                    
                    map_type = 0
                    
		    # All right if TID different or orientation different. Won't make any difference in those situations.	
                    if (al1.reference_start <= al2.reference_start):
			OUTER_D = abs(al1.reference_start - al2.reference_end)
			left_tid = al1.reference_name
			right_tid = al2.reference_name
			l_orient = al1.is_reverse
			r_orient = al2.is_reverse
			lf_target = al1.query_alignment_end # target mapping of left read in reference
			rt_target = al2.query_alignment_start # of right read in reference
			# sequencing orientation is always FR
			# assumes reads do not swap in reference except in case below (TD). Inverted insertion breakpoint margins will be overestimated.

			# reads likely swapped in alignment to reference
			if l_orient and not r_orient:
				
				lf_target = al1.query_alignment_start
				rt_target = al2.query_alignment_end

		    else:
			OUTER_D = abs(al2.reference_start - al1.reference_end)
			left_tid = al2.reference_name
			right_tid = al1.reference_name
			l_orient = al2.is_reverse
			r_orient = al1.is_reverse
			lf_target = al2.query_alignment_end
                        rt_target = al1.query_alignment_start

			if l_orient and not r_orient:

                                lf_target = al2.query_alignment_start
                                rt_target = al1.query_alignment_end
      
		CLType = str(int(l_orient)) + str(int(r_orient))
               
                # If discordant (due to mapping distance or orientation or mapping TID)
                # Check just to be safe even though discordant flag set in creating BAM files
                if map_type == 1 or (map_type == 0 and abs(OUTER_D-MEAN_D) > 3*SIG_D) or (map_type==0 and CLType!="01") or (map_type == 0 and left_tid != right_tid):
                        
                        temp=Cluster()
                        
                        if map_type == 0:
                           
			    # "Left" TID by convention is the one with the lower mapped coordinate if TIDs are different
                            if (al1.reference_start <= al2.reference_start) and (not al1.is_reverse) and (not al2.is_reverse):
                        
                                lf_bound = al1.reference_end                
                                rt_bound = al2.reference_end
      
                            elif (al1.reference_start <= al2.reference_start) and (not al1.is_reverse) and (al2.is_reverse):

                                lf_bound = al1.reference_end                
                                rt_bound = al2.reference_start
              			
                            elif (al1.reference_start <= al2.reference_start) and (al1.is_reverse) and (not al2.is_reverse):

                                lf_bound = al1.reference_start               
                                rt_bound = al2.reference_end
     
                            elif (al1.reference_start <= al2.reference_start) and (al1.is_reverse) and (al2.is_reverse):

                                lf_bound = al1.reference_start                
                                rt_bound = al2.reference_start
             			
                            elif (al1.reference_start > al2.reference_start) and (not al1.is_reverse) and (not al2.is_reverse):
                                
                                lf_bound = al2.reference_end              
                                rt_bound = al1.reference_end

                            elif (al1.reference_start > al2.reference_start) and (not al1.is_reverse) and (al2.is_reverse):

                                lf_bound = al2.reference_start               
                                rt_bound = al1.reference_end
         			
                            elif (al1.reference_start > al2.reference_start) and (al1.is_reverse) and (not al2.is_reverse):

                                lf_bound = al2.reference_end                
                                rt_bound = al1.reference_start
                                
                            elif (al1.reference_start > al2.reference_start) and (al1.is_reverse) and (al2.is_reverse):

                                lf_bound = al2.reference_start                
                                rt_bound = al1.reference_start
            				
                            temp.l_bound = lf_bound
                            temp.r_bound = rt_bound
                            temp.l_bound_O = lf_bound
                            temp.r_bound_O = rt_bound
                            temp.ltid = left_tid
                            temp.rtid = right_tid
						

                        if map_type == 1:
			    temp.l_bound = indivPos
			    temp.r_bound = -1
                            temp.ltid = left_tid
			    CLType[0] =  indivOrient
			    CLType[1] = 2
			    	
                            
                        temp.C_type = CLType
                        temp.Cmap_type = map_type
			temp.l_target = lf_target
			temp.r_target = rt_target

                        if temp.l_bound == None and temp.r_bound == None and temp.indiv_pos == None:
                            print "WARNING STORING BAD DATA ", temp


                        if map_type==0:
                            	
				DList1.append(temp)
 
                        if map_type == 1:
                            DList2.append(temp)

                # If even 1 alignment is concordant, don't use any mappings from read
                # Being safe, as discordant flag is anyway set when creating the input BAM files
                else:
		    DList1 =[]
		    DList2 = []	
                    return
			       
    
def ReadNextReadAlignments(bamname):

    bamfile = pysam.AlignmentFile(bamname)

    alignments = []
    qname = None

    for alignment in bamfile:
        if qname == None:
            qname = alignment.qname
            alignments.append(alignment)
        elif qname == alignment.qname:
            alignments.append(alignment)
        else:
            yield qname,alignments
            alignments = [alignment]
            qname = alignment.qname

    if qname != None:
        yield qname,alignments
                        
    
if __name__ == "__main__":

    f1 = open("../results/text/All_Discords_P.txt","w")
    f2 = open("../results/text/All_Discords_I.txt","w")
    
    currentFrag = 1

    # Read discordant alignments and write all possible discordant pairs to file. 

    #fp = open("../results/text/DiscordFxnTimeLog.txt","w")
    
    print "Entering loop..."

    start_for = time.clock()
    
    for (q1,aln1s),(q2,aln2s) in izip(ReadNextReadAlignments(FILE1),ReadNextReadAlignments(FILE2)):

        assert q1[:-2] == q2[:-2]

        DList1 = []
        DList2 = []
        
        start_fd = time.clock()

        FormDiscordant(aln1s, aln2s, DList1, DList2)
      
        for item in DList1:
             f1.write("%s %s\n" %(currentFrag, item))
        for item in DList2:
             f2.write("%s %s\n" %(currentFrag, item))
        

        end_fd = time.clock()
    
        currentFrag = currentFrag + 1

        #fp.write("Loop %d %f\n" %(currentFrag,(end_fd-start_fd)))   

    end_for = time.clock()

    print "Done with loop"

    f1.close()
    f2.close()
    #fp.close()

    #f=open("../results/text/TimeLog.txt","w")
    #f.write("For loop time: %f, form-discordant fxn time: %f\n" %(end_for-start_for,end_fd-start_fd))
    #f.close()
